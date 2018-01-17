function [particle_weights, gamma_dif, weight_dif, log_lh_time, graph_time, eig_time, m, aggregate_error_ratio] = LADelaunayLikelihoodFlexibleM(x_predicted, F, D, obs)
%   Function to compute the approximate posterior particles weights
%   The log-likelihood is computed in a distributed manner using Laplacian
%   approximation methods
%   In this method, the number of retained eigenvectors is not fixed but
%   set dynamically based on a pre-determined error threshold
%
%   Inputs:
%       x_predicted: (d+1)-by-N matrix of particle states, last row
%       corresponds to particle weights
%       F: Struct containing filter parameters
%       D: Struct containing measurement data
%       obs: Struct containing measurement model paraleters
%
% Output:
%       particle_weights: 1-by-N row vector of particle log-likelihood
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 14th, 2017

d = size(x_predicted,1)-1;

% First have each sensor compute local log-likelihood using only local
% measurements
log_lh_tic = tic;
for i=1:numel(D.sensorID)   
    z_received = D.measurements(:,i);
    % Compute expected measurement
    z_expected = obs.model(x_predicted(1:d,:), D.sensorLoc(:,i), obs);
    
    % Compute the Gaussian log-likelihood
    z_dif = F.minus(z_received, z_expected);
    
    log_lh_ss(i,:) = log(mvnpdf(z_dif', obs.mu', obs.R)+realmin)';
end
log_lh_time = toc(log_lh_tic);

% Construct the Delaunchy triangulation graph for all the particles
graph_tic = tic;
A = DelaunayGraph(x_predicted(1:2,:)');
% Construct Laplacian matrix
L = diag(sum(A,2)) - A;
graph_time = toc(graph_tic);

% Do eigenvalue decomposition of Laplacian matrix
eig_time_tic = tic;
[V_full,~] = eig(L);
eig_time = toc(eig_time_tic);

% Select the value for m
for m=1:F.LA.max_m
    log_m = (V_full(:,1:m)*V_full(:,1:m)'*(log_lh_ss'));
    log_N = (V_full*V_full'*(log_lh_ss'));
    
    lh_m = exp(log_m);
    lh_m = lh_m - max(lh_m);
    lh_m = lh_m./sum(lh_m);
    
    lh_N = exp(log_N);
    lh_N = lh_N - max(lh_N);
    lh_N = lh_N./sum(lh_N);
    
    lh_dif = sqrt(sum((lh_m - lh_N).^2,1));
    
    if sum(lh_dif>F.LA.max_lh_dif)==0
        break;
    end
end

% Select the m smallest eigenvectors;
V = V_full(:,1:m);

% Compute local coefficients
alpha_ss = V'*log_lh_ss';

% Sum up the local coefficients
if (F.gossip)
    [alpha, aggregate_error_ratio] = computeAggregateGossip(alpha_ss', F.A, F.max_gossip_iter);
    alpha = alpha';
else
    alpha = sum(alpha_ss,2);
    aggregate_error_ratio = zeros(1,F.LA.m);
end

% Compute approximate global joint log-likelihood
gamma_approx = (V*alpha)';

gamma = gamma_approx-max(gamma_approx);

% Compute unnormalized posterior weight
particle_weights = exp(gamma).*x_predicted(d+1,:);

% Give all particles equal weights if all particles have zero weight
if (sum(particle_weights) == 0)
    % This should never happen
    warning('All particle weights vanished');
    particle_weights = ones(1,N);
end

% Normalize the weights
particle_weights = particle_weights./sum(particle_weights); 

% Debug part
gamma_exact = sum(log_lh_ss);

gamma_dif = gamma_approx-gamma_exact;

gamma_exact = gamma_exact-max(gamma_exact);
weight_exact = exp(gamma_exact).*x_predicted(d+1,:);
weight_exact = weight_exact/sum(weight_exact);

weight_dif = weight_exact-particle_weights;