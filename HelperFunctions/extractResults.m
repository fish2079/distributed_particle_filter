% function [RMSE_trial, RMSE_time, RMSEFull, runtime, steptime, steptimeFull, N_effFull, detail] = extractResults(results, parameters)
function [RMSEFull, timeFull, weight_dif_full, N_effFull, aggregate_error_ratio, details] = extractResults(results, parameters)
%   Function to extract the simulation results
%
%   Input:
%       results: struct containing the tracking results
%       parameters: struct containing the parameters used in the simulatins

% initialize output variable
timeFull = [];
details = struct();


% % Average RMSE for each trial
% % The RMSE is stored as Nb_alg x Nb_trials matrix 
% RMSE_trial = reshape(mean(results.pos_error,2),[size(results.pos_error,1),size(results.pos_error,3)]);
% % Average RMSE over time
% % The RMSE is stored as Nb_alg x Nb_time_steps matrix
% RMSE_time = mean(results.pos_error,3);

% RMSE for each time step and each trial
RMSEFull = results.pos_error;

% Total runtime for each algorithm and each trial
% Nb_alg x Nb_trials matrix
% runtime = results.runtime;

% Extra results from the details subfield
% Loop through each trial
for tr=1:parameters.no_trials
    % Loop through each algorithm
    for alg = 1:numel(parameters.algorithms)
        try
            timeFull(alg,:,tr) = results.details{tr}{alg}.step_time;
            N_effFull(alg,:,tr) = results.details{tr}{alg}.N_eff;
            if (~isequal(func2str(parameters.algorithms{alg}),'BSpf'))
                weight_dif_full(alg,:,tr,:) = results.details{tr}{alg}.weight_dif;
                if (parameters.F.gossip)
                    aggregate_error_ratio(alg, :, tr) = mean(results.details{tr}{alg}.aggregate_error_ratio,2);
                else
                    aggregate_error_ratio(alg, :, tr) = zeros(1,50);
                end
            end
            switch func2str(parameters.algorithms{alg})
%                 case 'CSSpf_distributed'
%                     detail.CSSpf.weight_dif(alg, :, :, tr) = results.details{tr}{alg}.weight_dif;
                case 'LCpf_distributed'
                    details.LCpf.Hx_ss_dif(:,tr,:,:) = permute(results.details{tr}{alg}.Hx_ss_dif,[3,2,1]);
                case 'LApf_distributed'
                    details.LApf.gamma_dif(:,tr,:) = results.details{tr}{alg}.gamma_dif;
                    details.LApf.log_lh_time(:,tr) = results.details{tr}{alg}.log_lh_time;
                    details.LApf.graph_time(:,tr) = results.details{tr}{alg}.graph_time;
                    details.LApf.eig_time(:,tr) = results.details{tr}{alg}.eig_time;
                case 'LADelaunaypf_distributed'
                    details.LADelaunaypf.gamma_dif(:,tr,:) = results.details{tr}{alg}.gamma_dif;
                    details.LADelaunaypf.log_lh_time(:,tr) = results.details{tr}{alg}.log_lh_time;
                    details.LADelaunaypf.graph_time(:,tr) = results.details{tr}{alg}.graph_time;
                    details.LADelaunaypf.eig_time(:,tr) = results.details{tr}{alg}.eig_time;
                case 'Clusterpf_distributed'
                    details.Clusterpf.log_lh_time(:,tr) = results.details{tr}{alg}.log_lh_time';
                    details.Clusterpf.graph_time(:,tr) = results.details{tr}{alg}.graph_time';
                    details.Clusterpf.cluster_time(:,tr) = results.details{tr}{alg}.cluster_time';
                    details.Clusterpf.gamma_time(:,tr) = results.details{tr}{alg}.gamma_time';
                case 'ClusterDelaunaypf_distributed'
                    details.ClusterDelaunaypf.log_lh_time(:,tr) = results.details{tr}{alg}.log_lh_time';
                    details.ClusterDelaunaypf.graph_time(:,tr) = results.details{tr}{alg}.graph_time';
                    details.ClusterDelaunaypf.cluster_time(:,tr) = results.details{tr}{alg}.cluster_time';
                    details.ClusterDelaunaypf.gamma_time(:,tr) = results.details{tr}{alg}.gamma_time';
                case 'Debugpf'
                    details.Debugpf.weight_error(:,:,tr) = results.details{tr}{alg}.weight_error;
            end
%             steptime(alg,:) = mean(steptimeFull,3);
        end
    end
end

yo=5;
