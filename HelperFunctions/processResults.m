function [RMSE, time] = processResults(filename, mode)
%   Function to compute the average RMSE of position estimate and run time
%
%   Input:
%       filename: a cell array of strings, each cell contains one filename
%       Each file contains the tracking results for one specific simulation
%       setting (i.e., different values of N)
%       mode: flag to compute the average over time (1) or over trials (2)

% initialize output variable
averageError = [];
averageTime = [];

% Loop through all results
for i=1:numel(filename)
    load(filename{i});
    if (mode==1)
        % Average RMSE over time
        % The RMSE is stored as Nb_alg x Nb_trials matrix 
        % We concatenante the results from all files
        averageError = [averageError; reshape(mean(results.pos_error,2),[size(results.pos_error,1),size(results.pos_error,3)])];
    else
        % Average RMSE over trials
        % The RMSE is stored as Nb_alg x Nb_time matrix 
        % We concatenante the results from all files
        averageError = [averageError; reshape(mean(results.pos_error,3),[size(results.pos_error,1),size(results.pos_error,2)])];
    end
    % Compute average runtime
    % The runtime is stored as Nb_alg x Nb_trials matrix 
    % We concatenante the results from all files 
    averageTime = [averageTime; results.runtime];   
end

RMSE = averageError;
time = averageTime;