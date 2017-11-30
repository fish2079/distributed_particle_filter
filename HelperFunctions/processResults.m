function [RMSE, runtime, detail_time] = processResults(filename, mode)
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

if(nargin==1)
    mode = 1;
end

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
    
    % Compute total runtime for each trial
    runtime = results.runtime;
    
    
%     time = zeros(1, parameters.no_trials);
%     for j=1:parameters.no_trials
%         if(isfield(results.details{j}{1},'cluster_time'))
%             time(j) = sum(results.details{j}{1}.cluster_time)+sum(results.details{j}{1}.KNN_time)+sum(results.details{j}{1}.gamma_time);
%         else
%             time(j) = sum(results.details{j}{1}.step_time);
%         end
%     end
%     averageTime = [averageTime; time];   
end

RMSE = averageError;
runtime = averageTime;