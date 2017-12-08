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
sigma_vector = 5; %[1, 2.5, 5, 7.5, 10];

% Number of random trials
sim_parameters.no_trials = 100; 


% Flag for parallel run
sim_parameters.parallel = true;

% Flag for visualizing at each time step
sim_parameters.visualizeParticles = false;

% Tracking algorithms are
% 1. centralized bootstrap PF: BS
% 2. distributed CSS PF 
% 3. distributed LC PF
% 4. distributed Graph PF
alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LApf_distributed, @Clusterpf_distributed};
sim_parameters.algorithms = alg_lists([4]);

% Loop through each choice of particle number
for i=1:numel(sigma_vector)
    % Set number of particles
    sim_parameters.sigma = sigma_vector(i); 
    
    % Run the simulated track with all selected tracking algorithms 
    % Each filter uses N particles   
    [results, parameters]= runSimulatedTrack(sim_parameters);

    % Store the tracking results
    filename{i} = ['Track3_sigma',num2str(sim_parameters.sigma*10),'_trials',num2str(sim_parameters.no_trials),'.mat'];
    save(filename{i}, 'results','parameters');
end

% Plot the results
% plotRMSE(filename);