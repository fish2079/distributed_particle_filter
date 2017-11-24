function [RMSE, time] = processResults(filename)
%   Function to compute the average RMSE of position estimate and run time
%
%   Input:
%       filename: a cell array of strings, each cell contains one filename
%       Each file contains the tracking results for one specific simulation
%       setting (i.e., different values of N)

% initialize output variable
averageError = [];
averageTime = [];

% Loop through all results
for i=1:numel(filename)
    load(filename{i});
    % Compute RMSE over time
    % The RMSE is stored as Nb_alg x Nb_trials matrix 
    % We concatenante the results from all files
    averageError = [averageError; reshape(mean(results.pos_error,2),[size(results.pos_error,1),size(results.pos_error,3)])];
    
    % Compute average runtime
    % The runtime is stored as Nb_alg x Nb_trials matrix 
    % We concatenante the results from all files 
    averageTime = [averageTime; results.runtime];   
end

RMSE = averageError;
time = averageTime;