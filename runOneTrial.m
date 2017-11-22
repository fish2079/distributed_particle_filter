function [x_t, pos_error, runtime] = runOneTrial(S, F, dynamic, obs, algorithms, trial_idx)
%   Function to run one single Monte Carlo trial
%   The function generates the measurements for the trial and runs all
%   specified tracking algorithms
%
%   Input:
%       S: struct containing all track-related parameters
%       F: struct containing all filter-related parameters
%       dynamic: struct of transition model
%       obs: struct of measurement model
%       algorithms: a list of function handels for the tracking filters
%       trial_idx: initial seed for RNG, necessary to enable reproduction
%       of results
%
%   Output:
%       x_t: d x nb_steps matrix of estimated target states where d is
%       the target state dimension
%       pos_error: 1 x nb_steps row vector of position estimation error
%       runtime: 1 x nb_algs row vector of total runtime for each algorithm 

% Generate measurements for the given trial
rng(trial_idx+1000);
D = generateSimulatedMeasurements(S, obs);    
F_trial = F;

x_t = zeros(F.d, S.nb_steps, numel(algorithms));
pos_error = zeros(numel(algorithms), S.nb_steps);
runtime = zeros(1, numel(algorithms));

% Loop through all tracking filters
for alg = 1:numel(algorithms)
    disp(['trial ',num2str(trial_idx), ' ', func2str(algorithms{alg})]);
    % Set random seed to enable reproduction of results
    rng(trial_idx);

    % Set the tracking filter
    F_trial.algorithm = algorithms{alg};

    % run the algorithm and store the results for the trial
    tic;
    x_t(:,:,alg) = runFilter(S, F_trial, D, dynamic, obs);
    runtime(alg) = toc;
    pos_error(alg,:) = computePositionError(x_t(:,:,alg), S.x_t);
end