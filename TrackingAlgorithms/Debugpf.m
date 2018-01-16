function [x_updated, details] = Debugpf(x_old, F, D, dynamic, obs, details)
%   This function implements one time step of the centralized bootstrap 
%   particle filter
%
%   Inputs:
%       x_old: (d+1)*N matrix of particles from previous time step where N
%              is the number of particles and d is the dimension of the  
%              state. The d+1 row is the particle weights
%       F: Struct containing filter parameters
%       D: Struct containing measurement data
%       dynamics: Struct containing dynamic model parameters
%       obs: Struct containing measurement model paraleters
%
%   Outputs:
%       x_updated: (d+1)*N matrix of updated particles where N is the number of
%                  particles and d is the dimension of the state. The d+1 
%                  row is the updated particle weights
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017
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
    particle_weights = GaussianLikelihood([x_predicted; x_old(d+1,:)], F, D, obs);
    
    [CSS_weights, CSS_AER] = CSSLikelihood([x_predicted; x_old(d+1,:)], F, D, obs);
    
    [LC_weights, LC_AER] = LCLikelihood([x_predicted; x_old(d+1,:)], F, D, obs);
    
    [LA_weights, ~, ~, ~, ~, ~, LA_AER] = LADelaunayLikelihood([x_predicted; x_old(d+1,:)], F, D, obs);
    
    [Cluster_weights, ~, ~, ~, ~, ~, ~, CSS_AER] = ClusterDelaunayLikelihood([x_predicted; x_old(d+1,:)], F, D, obs);
    
    error = [norm(particle_weights-CSS_weights), norm(particle_weights-LC_weights), norm(particle_weights-LA_weights), norm(particle_weights-Cluster_weights)]';
    AER = [nanmean(CSS_AER), nanmean(LC_AER), nanmean(LA_AER), nanmean(CSS_AER)]';
    
    if(isfield(details, 'weight_error'))
        details.weight_error = [details.weight_error, error];
    else
        details.weight_error = error;
    end
    
    if(isfield(details, 'AER'))
        details.AER = [details.AER, error];
    else
        details.AER = AER;
    end
    
    if (1/sum(particle_weights.^2)<F.N_eff)
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



