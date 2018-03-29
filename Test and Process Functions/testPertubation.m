% Test script for comparing the algorithms' performance in terms of number 
% of particles on the simulated track
% For each parameter value, the script generates the simulated track,
% measurements and runs the specified tracking algorithms for a number of
% MC trials
% Note that all algorithms can only track one single target
% Multiple targets with target birth/death, clutter measurements are not
% considered

% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

warning('off','all');
clear;clc;
% Set path to helper functions
addpath('./DynamicModels/');
addpath('./HelperFunctions/');
addpath('./MeasurementModels/');
addpath('./TrackingAlgorithms/');

% Number of particles for the filter
pertubation_vector = [0.01, 0.1, 1, 5, 10, 25];%, 50, 100]; %[1,5,10,25,50,100,200,500]; %[1, 5,10,25];%,[50,100,200,500];

% Number of random trials
sim_parameters.no_trials = 40; 

sim_parameters.max_gossip_iter = 100;

% Flag for parallel run
sim_parameters.parallel = true;

% Flag for visualizing at each time step
sim_parameters.visualizeParticles = false;

% Flag for using gossip or exact aggregate
sim_parameters.gossip = false;

% Select the track
sim_parameters.track = 3;

%%
sim_parameters.nbEig = 9;

sim_parameters.nbClusters = 9; 

sim_parameters.max_degree = 2;
%%

% Tracking algorithms are
% 1. centralized bootstrap PF: BS
% 2. distributed CSS PF 
% 3. distributed LC PF
% 4. distributed Graph PF
alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LADelaunaypf_distributed, @ClusterDelaunaypf_distributed, @Debugpf, @LADelaunayFlexibleMpf_distributed};
% alg_lists = {@Debugpf};
sim_parameters.algorithms = alg_lists([3]);

% Loop through each choice of particle number
for i=1:numel(pertubation_vector)
    % Set number of particles
    sim_parameters.perturbation = pertubation_vector(i)/100;
    
    % Run the simulated track with all selected tracking algorithms 
    % Each filter uses N particles   
    [results, parameters]= runSimulatedTrack(sim_parameters);

   % Store the tracking results
    filename{i} = ['Track',num2str(sim_parameters.track)]; 
    filename{i} = [filename{i}, '_perturbation',num2str(parameters.perturbation*100)];
    filename{i} = [filename{i},'_N',num2str(parameters.F.N)];
    filename{i} = [filename{i},'_trials',num2str(parameters.no_trials)];
    filename{i} = [filename{i},'.mat'];
    save(filename{i}, 'results','parameters');
    save(filename{i}, 'results','parameters');
end

weight_error = [];
AER = [];
Neff = [];
for tr=1:40
    weight_error = cat(3, weight_error, results.details{tr}{1}.weight_error);
    AER = cat(3, AER, results.details{tr}{1}.AER);
    Neff = cat(3, Neff, results.details{tr}{1}.Neff);
end
w_error= [w_error; mean(mean(weight_error,3),2)'];
A_error = [A_error; mean(mean(AER,3),2)'];
% mean(mean(Neff,3),2)

