function [results, parameters]= runSimulatedTrack(sim_parameters)
%   This is the main function for running tests on simulated track.
%   The outputs from individual Monte Carlo trials are returned as output
%
%   Input:
%       sim_parameters: struct that contains number of trials, tracking
%       algorithms and number of particles
%           no_trials: number of random Monte Carlo trials
%           N: number of particles
%           algorithm: integer indicating the tracking algorithm
%                  1. centralized bootstrap PF
%                  2. distributed CSS PF
%                  3. distributed LC PF
%                  4. distributed Graph PF
%
%   Output:
%       results: struct array containing the results for each individual
%       trial
%       paramters: struct containing the parameters used in the tracking
%       scenario and algorithm
%
%   This function is designed to run random trials of simulated tracks in
%   parallel
%   The function generates the track, sensors and corresponding noise-free
%   measurements
%   Each trial is run on parallel with different measurement noise and
%   different algorithms are applied to the same dataset in each trial

% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th 2017

% Set the RNG seed 
% We always use the same random seed for track and measurements generation
% for consistency
rng('default');

% The struct dynamic contains all required parameters of target transition
% model
dynamic.model = @propagate_cv_ct; % dynamic model for target propagation
dynamic.a = 0.5; % turning rate
dynamic.p = 0.95; % probability of turning
dynamic.T = 1; % sampling interval 
dynamic.sigma_a = 0.05; % process noise standard deviation

% The struct S contains all relevant parameters for the target track
S.nb_steps = 50; % total number of time steps
S.initial = [45,45,5,0]';
S.area_length = 75; % tracking area size = area_length^2 (km^2)
S.nb_sensors = 9; % number of sensors
S.grid_sensor = true; % boolean flag to put the sensors in a grid
% Generate the track, sensors and corresponding noise-free measurements
S = buildSimulatedTrack(S, dynamic);

% The struct obs contains all required parameters of measurement model
obs.model = @computeBearing; % measurement model
obs.mu = 0;% [0,0]'; % measurement noise mean
obs.R = (5/180*pi)^2; %diag([5/180*pi, 0.005].^2);%diag([5^2,0.005^2]); % measurement noise covariance matrix
obs.sensorPos = S.sensorPos;
% Doppler specific parameters
obs.C = 0.343; % speed of sound in air in km/s
obs.F = 100;
obs.lambda = obs.F/obs.C;

% The struct F contains all relevant
% no of particles
if(isfield(sim_parameters, 'N'))
    F.N = sim_parameters.N; 
else
    F.N = 1000;
end
F.d = 4; % state vector dimension
% Particle Filter Regularization
F.sigma_noise = 10^-6; % std of Gaussian noise to add to particles at each iteration
F.initial_mu = S.x_t(:,1)'; % mean of initial particle cloud
F.initial_R = diag([0.5 0.5 0.05 0.05].^2); % covariance matrix of initial cloud
% Function handle for initial particle cloud generation
F.initializeState = @GaussianParticleCloudInitialization;
% Function handle for estimating target state from particle cloud
F.mmseStateEstimate = @mmseEstimateFromParticles;
% Bootstrap specific parameters
F.minus = @(z_exp, z_true) [wrapToPi(z_exp(1,:)-z_true(1,:))];%z_exp(2,:)-z_true(2,:)]; 
% LC specific parameters
[basis_x,basis_y]=meshgrid([-20:20:120],[20:20:160]);
F.LC.basis = [basis_x(:),basis_y(:)]';
F.LC.R = diag([10,10].^2);
% LA specific parameters
F.LA.KNN = 5; % number of nearest neighbors
F.LA.m = 10; % number of Eigenvectors to retain
% Clustering specific parameters
F.cluster.k = 5;
F.cluster.KNN = 5;

% Store the three parameter structs in the output struct
% The struct S contains parameters related to target track
% The struct F contains parameters related to tracking algorithm
parameters.S = S;
parameters.F = F;
parameters.no_trials = sim_parameters.no_trials; 
parameters.algorithms = sim_parameters.algorithms;

% Run all trials in parallel
parfor i = 1:parameters.no_trials
    % Run one single trial
    [x_t(:,:,:,i), pos_error(:,:,i), runtime(:,i)] = runOneTrial(S, F, dynamic, obs, parameters.algorithms, i);
end

% Store the results
results.x_t = x_t;
results.pos_error = pos_error;
results.runtime = runtime;
