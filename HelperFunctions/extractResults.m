function [RMSE, RMSEFull, runtime, runtimeFull, detail] = extractResults(filename)
%   Function to compute the average RMSE of position estimate and run time
%
%   Input:
%       filename: a cell array of strings, each cell contains one filename
%       Each file contains the tracking results for one specific simulation
%       setting (i.e., different values of N)

% initialize output variable
runtimeFull = [];
detail = struct();

% Load the matfile
load(filename);

% Average RMSE over time
% The RMSE is stored as Nb_alg x Nb_trials matrix 
% We concatenante the results from all files
RMSE = reshape(mean(results.pos_error,2),[size(results.pos_error,1),size(results.pos_error,3)]);

% RMSE for each time step and each trial
RMSEFull = results.pos_error;

% Average runtime over time
runtime = results.runtime;

runtimeFull = [];
% Extra results from the details subfield
% Loop through each trial
for tr=1:parameters.no_trials
    % Loop through each algorithm
    for alg = 1:numel(parameters.algorithms)
        if (isfield(results, 'details'))
            runtimeFull(alg,:,tr) = results.details{tr}{alg}.step_time;
            try
                switch func2str(parameters.algorithms{alg})
                    case 'CSSpf_distributed'
                        detail.CSSpf.weight_error(:,tr) = results.details{tr}{alg}.weight_error';
                    case 'LCpf_distributed'
                        detail.LCpf.Hx_ss_dif(:,tr,:) = results.details{tr}{alg}.Hx_ss_dif;
                    case 'LApf_distributed'
                        detail.LApf.gamma_dif(:,tr) = results.details{tr}{alg}.gamma_dif;
                        detail.LApf.weight_dif(:,tr) = results.details{tr}{alg}.weight_dif;
                        detail.LApf.KNN_time(:,tr) = results.details{tr}{alg}.KNN_time;
                        detail.LApf.eig_time(:,tr) = results.details{tr}{alg}.eig_time;
                    case 'Clusterpf_distributed'
                        detail.Clusterpf.KNN_time(:,tr) = results.details{tr}{alg}.KNN_time';
                        detail.Clusterpf.cluster_time(:,tr) = results.details{tr}{alg}.cluster_time';
                        detail.Clusterpf.gamma_time(:,tr) = results.details{tr}{alg}.gamma_time';
                    case 'Debugpf'
                        detail.Debugpf.weight_error(:,:,tr) = results.details{tr}{alg}.weight_error;
                end
            end
        end
    end
end
