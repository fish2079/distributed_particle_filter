function testN(track)
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
% Set path to helper functions
addpath('./DynamicModels/');
addpath('./HelperFunctions/');
addpath('./MeasurementModels/');
addpath('./TrackingAlgorithms/');

N_vector = 100:100:2000;

% Number of random trials
sim_parameters.no_trials = 200; 

switch track
    case 2
        sim_parameters.measModel = 'bearing';
    case 3
        sim_parameters.measModel = 'range';
end

% Flag for parallel run
sim_parameters.parallel = true;

% Flag for visualizing at each time step
sim_parameters.visualizeParticles = false;

% Flag for using gossip or exact aggregate
sim_parameters.gossip = false;

% sim_parameters.graphMethod = 'Delaunay';
% 
% sim_parameters.weightedEdge = 1;
% sim_parameters.weightedEdgeStyle = 1;

% Tracking algorithms are
% 1. centralized bootstrap PF: BS
% 2. distributed CSS PF 
% 3. distributed LC PF
% 4. distributed Graph PF
alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LADelaunaypf_distributed, @ClusterDelaunaypf_distributed};
sim_parameters.algorithms = alg_lists(1);

% Select the track
sim_parameters.track = track;
for i=1:numel(N_vector)
    sim_parameters.N = N_vector(i);
    [results, parameters]= runSimulatedTrack(sim_parameters);
    
    filename{i} = ['Track',num2str(track),'_BSpf'];
    filename{i} = [filename{i},'_N',num2str(parameters.F.N)];
    filename{i} = [filename{i},'_trials',num2str(parameters.no_trials)];
    filename{i} = [filename{i},'.mat'];
    save(filename{i}, 'results','parameters');
end