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
% dynamic.a = 0.5; % turning rate
dynamic.a = -0.25; % turning rate track 2

dynamic.p = 0.05; % probability of turning tracke 1
dynamic.p = 0.75; % probability of turning tracke 2
% dynamic.p = 0.95; % probability of turning track 3
dynamic.T = 1; % sampling interval 
dynamic.sigma_a = 0.05; % process noise standard deviation

% The struct S contains all relevant parameters for the target track
S.nb_steps = 50; % total number of time steps
S.initial = [45,45,5,0]';
if(isfield(sim_parameters,'areaLength'))
    S.area_length = sim_parameters.areaLength;
else
    S.area_length = 75; % tracking area size = area_length^2 (km^2)
end

% S.nb_sensors = 4; % number of sensors % Track 1
S.nb_sensors = 16; % number of sensors % Track 2
% S.nb_sensors = 9; % number of sensors Track 3

S.grid_sensor = true; % boolean flag to put the sensors in a grid
% Generate the track, sensors and corresponding noise-free measurements
S = buildSimulatedTrack(S, dynamic);
S.parallel = sim_parameters.parallel;
S.visualizeParticles = sim_parameters.visualizeParticles;

% The struct obs contains all required parameters of measurement model
obs.model = @computeBearing; % measurement model
obs.mu = 0;% [0,0]'; % measurement noise mean
if (isfield(sim_parameters,'sigma'))
    obs.R = (sim_parameters.sigma/180*pi)^2;
else
    obs.R = (5/180*pi)^2; %diag([5/180*pi, 0.005].^2);%diag([5^2,0.005^2]); % measurement noise covariance matrix
end
obs.sensorPos = S.sensorPos;
% Doppler specific parameters
% obs.C = 0.343; % speed of sound in air in km/s
% obs.F = 100;
% obs.lambda = obs.F/obs.C;

% The struct F contains all relevant
% no of particles
if(isfield(sim_parameters, 'N'))
    F.N = sim_parameters.N; 
else
    F.N = 500;
end
F.N_eff = F.N/3; % minimum number of effective particles before resampling
F.d = 4; % state vector dimension
if (isfield(sim_parameters,'gossip'))
    F.gossip = sim_parameters.gossip;
else
    F.gossip = false;
end
F.A = S.A; % adjacency matrix of the sensor network

if (isfield(sim_parameters,'max_gossip_iter'))
    F.max_gossip_iter = sim_parameters.max_gossip_iter;
else
    F.max_gossip_iter = 50; % maximum number of gossip iterations
end

% Particle Filter Regularization
F.sigma_noise = 10^-3; % std of Gaussian noise to add to particles at each iteration
F.initial_mu = S.x_t(:,1)'; % mean of initial particle cloud
F.initial_R = diag([0.5 0.5 0.05 0.05].^2); % covariance matrix of initial cloud
% Function handle for initial particle cloud generation
F.initializeState = @GaussianParticleCloudInitialization;
% Function handle for estimating target state from particle cloud
F.mmseStateEstimate = @mmseEstimateFromParticles;
% Bootstrap specific parameters
F.minus = @(z_exp, z_true) [wrapToPi(z_exp(1,:)-z_true(1,:))];%z_exp(2,:)-z_true(2,:)]; 

% LC specific parameters
% [basis_x,basis_y]=meshgrid([-20:20:120],[20:20:160]);
% F.LC.basis = [basis_x(:),basis_y(:)]';
% F.LC.R = diag([10,10].^2);
if(isfield(sim_parameters, 'max_degree'))
    F.LC.max_degree = sim_parameters.max_degree;
else   
    F.LC.max_degree = 1;
end

% LA specific parameters
% Number of neighbors is irrelevant since we now use Delaunay
% triangulation
% Maintain the parameters just in case
if(isfield(sim_parameters, 'KNN'))
    F.LA.KNN = sim_parameters.KNN; % number of nearest neighbors
else
    F.LA.KNN = 10; 
end

if(isfield(sim_parameters, 'nbEig'))
    F.LA.m = sim_parameters.nbEig; % number of Eigenvectors to retain
else
    F.LA.m = 6;
end

% Clustering specific parameters
if(isfield(sim_parameters,'nbClusters'))
    F.cluster.k = sim_parameters.nbClusters;
else
    F.cluster.k = 6;
end

% As in the case LApf, KNN is obsolete since we use Delaunay triangulation
% only
F.cluster.KNN = F.LA.KNN;

% Store the two parameter structs in the output struct
% The struct S contains parameters related to target track
% The struct F contains parameters related to tracking algorithm
parameters = sim_parameters;
parameters.S = S;
parameters.F = F;

% Run all trials in parallel if required
if (parameters.parallel)
    parfor i = 1:parameters.no_trials
        % Run one single trial
        [x_t(:,:,:,i), pos_error(:,:,i), runtime(:,i), details{i}] = runOneTrial(S, F, dynamic, obs, parameters.algorithms, i);
    end
else
    for i = 1:parameters.no_trials
        % Run one single trial
        [x_t(:,:,:,i), pos_error(:,:,i), runtime(:,i), details{i}] = runOneTrial(S, F, dynamic, obs, parameters.algorithms, i);
    end
end

% Store the results
results.x_t = x_t;
results.pos_error = pos_error;
results.runtime = runtime;
results.details = details;