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

% The struct S contains all relevant parameters for the target track
[S, dynamic] = buildTrackParameters(sim_parameters);

% Generate the track, sensors and corresponding noise-free measurements
S = buildSimulatedTrack(S, dynamic);
S.parallel = sim_parameters.parallel;
S.visualizeParticles = sim_parameters.visualizeParticles;

% The struct obs contains all required parameters of measurement model
if (isfield(sim_parameters,'measModel'))
    obs.measModel = sim_parameters.measModel;
else
    obs.measModel = 'bearing';
end
    
switch obs.measModel
    case 'bearing'
        obs.model = @computeBearing; % measurement model
        obs.mu = 0; % measurement noise mean
        if (isfield(sim_parameters,'sigma'))
            obs.R = (sim_parameters.sigma/180*pi)^2;
        else
            obs.R = (5/180*pi)^2; % measurement noise covariance matrix
        end
        F.minus = @(z_exp, z_true) [wrapToPi(z_exp(1,:)-z_true(1,:))]; 
    case 'range'
        obs.model = @computeRange; % measurement model
        obs.mu = 0; % measurement noise mean
        if (isfield(sim_parameters,'sigma'))
            obs.R = sim_parameters.sigma^2;
        else
            obs.R = 5^2; % measurement noise covariance matrix
        end
        F.minus = @(z_exp, z_true) z_exp-z_true; 
    case 'bearDoppler'
        obs.model = @computeBearingDoppler;
        obs.mu = [0,0]';
        obs.C = 0.343; % speed of sound in air in km/s
        obs.F = 100;
        obs.lambda = obs.F/obs.C;
        if (isfield(sim_parameters,'sigma'))
            obs.R = diag(sim_parameters.sigma.^2);
        else
            obs.R = diag([5/180*pi,0.5].^2); % measurement noise covariance matrix
        end
        F.minus = @(z_exp, z_true) [wrapToPi(z_exp(1,:)-z_true(1,:)); z_exp(2,:)-z_true(2,:)]; 
end

obs.sensorPos = S.sensorPos;

% The struct F contains all relevant
% no of particles
if(isfield(sim_parameters, 'N'))
    F.N = sim_parameters.N; 
else
    F.N = 1000;
end
F.N_eff = F.N/3; % minimum number of effective particles before resampling
F.d = 4; % state vector dimension
if (isfield(sim_parameters,'gossip'))
    F.gossip = sim_parameters.gossip;
else
    F.gossip = false;
end
F.A = S.A; % adjacency matrix of the sensor network

if (isfield(sim_parameters,'perturbation'))
    F.perturbation = sim_parameters.perturbation;
else
    F.perturbation = 0;
end

if (isfield(sim_parameters,'max_gossip_iter'))
    F.max_gossip_iter = sim_parameters.max_gossip_iter;
else
    F.max_gossip_iter = 50; % maximum number of gossip iterations
end

% Particle Filter Regularization
F.sigma_noise = 10^-3; % std of Gaussian noise to add to particles at each iteration

if (isfield(sim_parameters,'initialization'))
    switch sim_parameters.initialization
        case 'perfect'
            F.initial_mu = S.x_t(:,1)'; % mean of initial particle cloud
            F.initial_R = diag([0.5 0.5 0.05 0.05].^2); % covariance matrix of initial cloud
            % Function handle for initial particle cloud generation
            F.initializeState = @GaussianParticleCloudInitialization;
        otherwise
            F.initial_mu = S.x_t(:,1)'; % mean of initial particle cloud
            F.initial_R = diag([0.5 0.5 0.05 0.05].^2); % covariance matrix of initial cloud
            F.initial_particles_R = diag([5 5 0.5 0.5].^2);
            % Function handle for initial particle cloud generation
            F.initializeState = @InaccurateParticleCloudInitialization;
    end
else
    F.initial_mu = S.x_t(:,1)'; % mean of initial particle cloud
    F.initial_R = diag([0.5 0.5 0.05 0.05].^2); % covariance matrix of initial cloud
    F.initial_particles_R = diag([5 5 0.5 0.5].^2);
    % Function handle for initial particle cloud generation
    F.initializeState = @InaccurateParticleCloudInitialization;
end
% F.initial_mu = S.x_t(:,1)'; % mean of initial particle cloud
% F.initial_R = diag([0.5 0.5 0.05 0.05].^2); % covariance matrix of initial cloud
% % Function handle for initial particle cloud generation
% F.initializeState = @GaussianParticleCloudInitialization;



% F.initial_lower_limit = S.x_t(:,1)'-[10,10,5,5]; 
% F.initial_upper_limit = S.x_t(:,1)'+[10,10,5,5]; 
% F.initializeState = @UniformParticleCloudInitialization;

% F.initial_sensor = S.sensorPos; 
% F.initializeState = @BearingCrossParticleCloudInitialization;

% Function handle for estimating target state from particle cloud
F.mmseStateEstimate = @mmseEstimateFromParticles;

% LC specific parameters
% [basis_x,basis_y]=meshgrid([-20:20:120],[20:20:160]);
% F.LC.basis = [basis_x(:),basis_y(:)]';
% F.LC.R = diag([10,10].^2);
if(isfield(sim_parameters, 'max_degree'))
    F.LC.max_degree = sim_parameters.max_degree;
else   
    F.LC.max_degree = 2;
end

if(isfield(sim_parameters, 'max_degree_GS'))
    F.LC.max_degree_GS = sim_parameters.max_degree_GS;
else   
    F.LC.max_degree_GS = F.LC.max_degree;
end

% LA specific parameters
% Graph construction methods
if(isfield(sim_parameters, 'graphMethod'))
    F.LA.graphMethod = sim_parameters.graphMethod;
    switch sim_parameters.graphMethod
        case 'KNN'
            % Number of neighbors for KNN graph
            if(isfield(sim_parameters, 'KNN'))
                F.LA.KNN = sim_parameters.KNN; % number of nearest neighbors
            else
                F.LA.KNN = 10; 
            end
        case 'Epsilon'
            % Number of neighbors for KNN graph
            if(isfield(sim_parameters, 'Epsilon'))
                F.LA.epsilon = sim_parameters.Epsilon; % number of nearest neighbors
            else
                F.LA.epsilon = 1/3; 
            end
    end
end

if(isfield(sim_parameters, 'weightedEdge'))
    F.LA.weightedEdge = sim_parameters.weightedEdge; % flag to use the weighted adjacency matrix
    F.LA.weightedEdgeStyle = sim_parameters.weightedEdgeStyle; % choose different weighted matrix
else
    F.LA.weightedEdge = false;
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

% Graph construction methods
if(isfield(sim_parameters, 'graphMethod'))
    F.cluster.graphMethod = sim_parameters.graphMethod;
    switch sim_parameters.graphMethod
        case 'KNN'
            % Number of neighbors for KNN graph
            if(isfield(sim_parameters, 'KNN'))
                F.cluster.KNN = sim_parameters.KNN; % number of nearest neighbors
            else
                F.cluster.KNN = 10; 
            end
        case 'Epsilon'
            % Number of neighbors for KNN graph
            if(isfield(sim_parameters, 'epsilon'))
                F.cluster.epsilon = sim_parameters.epsilon; % number of nearest neighbors
            else
                F.cluster.epsilon = 1/3; 
            end
    end
end

if(isfield(sim_parameters, 'weightedEdge'))
    F.cluster.weightedEdge = sim_parameters.weightedEdge; % flag to use the weighted adjacency matrix
    F.cluster.weightedEdgeStyle = sim_parameters.weightedEdgeStyle; % choose different weighted matrix
else
    F.cluster.weightedEdge = false;
end

% Store the three parameter structs in the output struct
% The struct S contains parameters related to target track
% The struct F contains parameters related to tracking algorithm
% The struct obs contains parameters related to observation model
parameters = sim_parameters;
parameters.S = S;
parameters.F = F;
parameters.obs = obs;

% Run all trials in parallel if required
if (parameters.parallel)
    parfor i = 1:parameters.no_trials
        % Run one single trial
        try
            [x_t(:,:,:,i), pos_error(:,:,i), runtime(:,i), details{i}] = runOneTrial(S, F, dynamic, obs, parameters.algorithms, i);
        catch ME
            i
            ME
        end
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