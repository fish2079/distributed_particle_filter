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
N_vector = 500; %[250, 500, 750, 1000];

% Number of random trials
sim_parameters.no_trials = 40; 

sim_parameters.max_gossip_iter = 50;

% Flag for parallel run
sim_parameters.parallel = true;

% Flag for visualizing at each time step
sim_parameters.visualizeParticles = false;

% Flag for using gossip or exact aggregate
sim_parameters.gossip = true;

% Select the track
sim_parameters.track = 2;

% Tracking algorithms are
% 1. centralized bootstrap PF: BS
% 2. distributed CSS PF 
% 3. distributed LC PF
% 4. distributed Graph PF
% alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LApf_distributed, @LADelaunaypf_distributed, @Clusterpf_distributed, @ClusterDelaunaypf_distributed};
alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LADelaunaypf_distributed, @ClusterDelaunaypf_distributed, @Debugpf};
sim_parameters.algorithms = alg_lists([1:5]);

% sim_parameters.areaLength = 125;
% Loop through each choice of particle number
for i=1:numel(N_vector)
    % Set number of particles
    sim_parameters.N = N_vector(i); 
    
    % Run the simulated track with all selected tracking algorithms 
    % Each filter uses N particles   
    [results, parameters]= runSimulatedTrack(sim_parameters);

    % Store the tracking results
    filename{i} = ['Track',num2str(sim_parameters.track)]; 
    filename{i} = [filename{i}, '_gossip',num2str(parameters.max_gossip_iter)];
    filename{i} = [filename{i},'_N',num2str(parameters.F.N)];
    filename{i} = [filename{i},'_trials',num2str(parameters.no_trials)];
    filename{i} = [filename{i},'.mat'];
    save(filename{i}, 'results','parameters');
    save(filename{i}, 'results','parameters');
end

mean(mean(results.pos_error,3),2)
mean(results.runtime,2)
weight_error = [];
AER = [];
for tr=1:20
    weight_error = cat(3, weight_error, results.details{tr}{1}.weight_error);
    AER = cat(3, AER, results.details{tr}{1}.AER);
end
mean(mean(weight_error,3),2)
mean(mean(AER,3),2)