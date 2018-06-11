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
maxDegree_vector = 2; %[1,3,4,5];
gossip_vector = [1, 5, 10, 25, 50, 75, 100, 150, 200];

% Number of random trials
sim_parameters.no_trials = 200; 

% Flag for parallel run
sim_parameters.parallel = true;

% Flag for visualizing at each time step
sim_parameters.visualizeParticles = false;

% Flag for using gossip or exact aggregate
sim_parameters.gossip = true;

% Select the track
sim_parameters.track = 3;
if (sim_parameters.track==2)
    sim_parameters.measModel = 'bearing';
else
    sim_parameters.measModel = 'range';
end

% Tracking algorithms are
% 1. centralized bootstrap PF: BS
% 2. distributed CSS PF 
% 3. distributed LC PF
% 4. distributed Graph PF
alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LCpf_GS_distributed, @LApf_distributed, @Clusterpf_distributed};
sim_parameters.algorithms = alg_lists(4);

% Loop through each choice of particle number
for i=1:numel(maxDegree_vector)
    % Set number of particles
    sim_parameters.max_degree = maxDegree_vector(i); 
    
    for j=1:numel(gossip_vector)
        sim_parameters.max_gossip_iter = gossip_vector(j);
        % Run the simulated track with all selected tracking algorithms 
        % Each filter uses N particles   
        [results, parameters]= runSimulatedTrack(sim_parameters);

        % Store the tracking results
        filename = ['Track',num2str(sim_parameters.track),'_LCpf-GS_'];
        filename = [filename, '_gossip',num2str(gossip_vector(j))];
        filename = [filename,'_maxDegree',num2str(parameters.F.LC.max_degree)];
        filename = [filename,'_N',num2str(parameters.F.N)];
        filename = [filename,'_trials',num2str(parameters.no_trials)];
        filename = [filename,'.mat'];
        save(filename, 'results','parameters');
    end
end

% Plot the results
% plotRMSE(filename);