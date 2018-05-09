warning('off','all');
clear;clc;
% Set path to helper functions
addpath('./DynamicModels/');
addpath('./HelperFunctions/');
addpath('./MeasurementModels/');
addpath('./TrackingAlgorithms/');

% Number of particles for the filter
N_vector = 500; 
k_vector = 10:10:500; 

% Tracking algorithms are
% 1. centralized bootstrap PF: BS
% 2. distributed CSS PF 
% 3. distributed LC PF
% 4. distributed Graph PF
alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LADelaunaypf_distributed, @ClusterDelaunaypf_distributed};
sim_parameters.algorithms = alg_lists(4);

% Select the track
sim_parameters.track = 3;

% Loop through each choice of particle number
for i=1:numel(k_vector)   
%     sim_parameters.KNN = KNN_vector(1);   
    
    % Generate target position
    x_t = [5,5]'; 
    
    obs.model = @computeRange; % measurement model
    obs.mu = 0; % measurement noise mean
    obs.R = 5^2; % measurement noise covariance matrix
    
    F.minus = @(z_exp, z_true) z_exp-z_true; 
    F.cluster.KNNgraph = false;
    F.cluster.KNN = 100;
    F.cluster.weightedEdge = false;
    F.N = 500;
    F.cluster.k = k_vector(i);
    F.gossip = false;
    F.perturbation = 0;
    
    % Generate target measurement
    range = computeRange(x_t, [0;0], obs)+normrnd(0,sqrt(obs.R));
    D.measurements = range;
    D.sensorLoc = [0,0]';
    D.sensorID = 1;
        
    parfor tr=1:200
        rng(tr);
        % Generate the particles
        % The particles are fixed across trials
        x_predicted = rand(2,500)*100-50;
        x_predicted(3:4,:) = 0;
        x_predicted(5,:) = 1/500;
    
        [~,~,weight_dif(tr)] = ClusterDelaunayLikelihood(x_predicted, F, D, obs);
    end
   
    weight_dif_matrix(i,:) = weight_dif;
end

% Store the tracking results
filename{i} = ['Track3_SingleTime_LApf'];
filename{i} = [filename{i},'_m',num2str(parameters.F.LA.m)];
filename{i} = [filename{i},'_N',num2str(parameters.F.N)];
filename{i} = [filename{i},'_trials',num2str(parameters.no_trials)];
filename{i} = [filename{i},'.mat'];
save(filename{i}, 'results','parameters');