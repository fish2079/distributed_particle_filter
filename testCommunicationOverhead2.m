% Test script for comparing the algorithms' performance in terms of communication
% overhead on the simulated track
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
N = 500; 
overhead_vector = [60:60:1200];

% Number of random trials
sim_parameters.no_trials = 200;

sim_parameters.max_gossip_iter = 30;

% Flag for parallel run
sim_parameters.parallel = true;

% Flag for visualizing at each time step
sim_parameters.visualizeParticles = false;

% Flag for using gossip or exact aggregate
sim_parameters.gossip = true;

% Tracking algorithms are
% 1. centralized bootstrap PF: BS
% 2. distributed CSS PF 
% 3. distributed LC PF
% 4. distributed Graph PF
% alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LApf_distributed, @LADelaunaypf_distributed, @Clusterpf_distributed, @ClusterDelaunaypf_distributed};
alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LADelaunaypf_distributed, @ClusterDelaunaypf_distributed};

% sim_parameters.algorithms = alg_lists([1:5]);

sim_parameters.N = N; 
sim_parameters.max_gossip_iter = sim_parameters.max_gossip_iter; 
for i=1:numel(overhead_vector)
    for alg=4:5
        sim_parameters.algorithms = alg_lists(alg);
        
        sim_parameters.nbEig = overhead_vector(i)/sim_parameters.max_gossip_iter;
        sim_parameters.nbClusters = overhead_vector(i)/sim_parameters.max_gossip_iter;
        

        % Run the simulated track with all selected tracking algorithms 
        % Each filter uses N particles   
        [results, parameters]= runSimulatedTrack(sim_parameters);

        % Store the tracking results
        filename{i} = ['Track1_'];
        filename{i} = [filename{i}, func2str(alg_lists{alg})];
        filename{i} = [filename{i}, '_overhead',num2str(overhead_vector(i))];
        filename{i} = [filename{i}, '_m',num2str(sim_parameters.nbEig)];
        filename{i} = [filename{i}, '_C',num2str(sim_parameters.nbClusters)];
        filename{i} = [filename{i},'_N',num2str(parameters.F.N)];
        filename{i} = [filename{i},'_trials',num2str(parameters.no_trials)];
        filename{i} = [filename{i},'.mat'];
        save(filename{i}, 'results','parameters');
    end
end