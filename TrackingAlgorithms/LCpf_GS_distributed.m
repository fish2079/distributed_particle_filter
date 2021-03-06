function [x_updated, details]= LCpf_GS_distributed(x_old, F, D, dynamic, obs, details)
%   This function implements one step of the distributed likelihood 
%   consensus particle filter. 
%
%   Inputs:
%       x_old: (d+1)*N matrix of particles from previous time step where N
%              is the number of particles and d is the dimension of the  
%              state. The d+1 row is the particle weights
%       F: Struct containing filter parameters
%       D: Struct containing measurement data
%       dynamics: Struct containing dynamic model parameters
%       obs: Struct containing measurement model paraleters
%       details: Struct containing debug data for the algorithm
%
%   Outputs:
%       x_updated: (d+1)*N matrix of updated particles where N is the number of
%                  particles and d is the dimension of the state. The d+1 
%                  row is the updated particle weights
%       details: updated input struct of the same name
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 13th, 2017

% Note: This code only maintains a single copy of the particles, rather
% than one for each node in the network. Because we assume/enforce that the
% weighted particle clouds remain precisely synchronized (by running max
% consensus after each distributed averaging update, and under the
% assumption of reliable communications) the weighted particle clouds at
% all nodes remain the same, so this approach is equivalent. It also uses
% less memory and is less computationally intensive.
step_tic = tic;
N = F.N; % number of particles
d = F.d; % state dimension

% Generate regularization noise
regularization_noise = F.sigma_noise*randn(d,N);

% Propagate particles if necessary
if(F.initial)
    x_predicted = x_old(1:d,:);
else
    x_predicted = dynamic.model(x_old(1:d,:), dynamic);
end

% Proceed if there are measurements 
if ~isempty(D.measurements) 
    % Compute the particle likelihood
%     [particle_weights, ~, aggregate_error_ratio] = LCLikelihoodMeasModel([x_predicted; x_old(d+1,:)], F, D, obs);
    [particle_weights, aggregate_error_ratio] = LCLikelihood_GS([x_predicted; x_old(d+1,:)], F, D, obs);
    Hx_ss_dif = 0;
   
    [bs_weights] = GaussianLikelihood([x_predicted; x_old(d+1,:)], F, D, obs);
    
    if (isfield(details,'weight_dif'))
        details.weight_dif = [details.weight_dif, norm(particle_weights-bs_weights)];
    else
        details.weight_dif = norm(particle_weights-bs_weights);
    end
    
    if (isfield(details,'Hx_ss_dif'))
        details.Hx_ss_dif = cat(3, details.Hx_ss_dif, Hx_ss_dif);
    else
        details.Hx_ss_dif = Hx_ss_dif;
    end
    
    if (F.gossip)
        if (isfield(details,'aggregate_error_ratio'))
            details.aggregate_error_ratio = [details.aggregate_error_ratio; aggregate_error_ratio];
        else
            details.aggregate_error_ratio = aggregate_error_ratio;
        end
    end
    
    N_eff = 1/sum(particle_weights.^2);
    if (N_eff<F.N_eff)
        % Sample according to weights with replacement
        I = randsample((1:N)', N, true, particle_weights);

        % Add regularization noise and set the weights
        x_updated = [ x_predicted(1:d,I) + regularization_noise; ones(1,N)/N ];
    else
        x_updated = [ x_predicted(1:d,:) + regularization_noise; particle_weights];
    end
else
    % If there is no measurement, propagate predicted particles and assign 
    % them equal weights
    x_updated = [ x_predicted + regularization_noise; ones(1,N)/N ];
end
if (isfield(details,'step_time'))
    details.step_time = [details.step_time, toc(step_tic)];
else
    details.step_time = toc(step_tic);
end

if (isfield(details,'N_eff'))
    details.N_eff = [details.N_eff, N_eff];
else
    details.N_eff = N_eff;
end
