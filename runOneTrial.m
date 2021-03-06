function [x_t, pos_error, runtime, details] = runOneTrial(S, F, dynamic, obs, algorithms, trial_idx)
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
%       details: 1 x nb_algs stuct array, each struct contains data that
%       can be used for debugging individual algorithms. 

% Generate measurements for the given trial
rng(trial_idx+1000);
D = generateSimulatedMeasurements(S, obs);    
F_trial = F;
F_trial.initial_sensors =  D.sensorLoc{1};
F_trial.initial_meas =  D.measurements{1};

x_t = zeros(F.d, S.nb_steps, numel(algorithms));
pos_error = zeros(numel(algorithms), S.nb_steps);
runtime = zeros(1, numel(algorithms));

% Loop through all tracking filters
for alg = 1:numel(algorithms)
    disp(['trial ',num2str(trial_idx), ' ', func2str(algorithms{alg})]);
    % Set random seed to enable reproduction of results
    rng(trial_idx, 'twister');

    % Set the tracking filter
    F_trial.algorithm = algorithms{alg};

    % run the algorithm and store the results for the trial
    trial_time=tic;
    [x_t(:,:,alg), details{alg}] = runFilter(S, F_trial, D, dynamic, obs);
    runtime(alg) = toc(trial_time);
    pos_error(alg,:) = computePositionError(x_t(:,:,alg), S.x_t);
end