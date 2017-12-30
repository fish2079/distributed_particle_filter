function [particle_weights, gamma_dif, weight_dif, cluster_time, log_lh_time, graph_time, gamma_time, aggregate_error_ratio] = ClusterDelaunayLikelihood(x_predicted, F, D, obs)
%   Function to compute the approximate posterior particles weights
%   The log-likelihood is computed in a distributed manner by clustering
%   particles and gossiping on cluster log-likelihoods
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
cluster_tic = tic;
% Start by clustering the particles into k clusters
[idx, ~] = kmeans(x_predicted',F.cluster.k);
cluster_time = toc(cluster_tic);

% Compute clustering assignment matrix
C = zeros(F.cluster.k, F.N);
ind = sub2ind([F.cluster.k,F.N],idx',1:F.N);
C(ind) = 1;

% Now have each sensor compute local log-likelihood using only local
% measurements
% Also, compute the log-likelihood of clusters by summing up the
% log-likelihoood of particles in that cluster
log_lh_tic = tic;
for i=1:numel(D.sensorID)   
    z_received = D.measurements(:,i);
    % Compute expected measurement
    z_expected = obs.model(x_predicted(1:d,:), D.sensorLoc(:,i), obs);
    
    % Compute the Gaussian log-likelihood
    z_dif = F.minus(z_received, z_expected);
    
    log_lh_ss(i,:) = log(mvnpdf(z_dif', obs.mu', obs.R))';
  
    log_lh_cluster_ss(i,:) = C*log_lh_ss(i,:)';
end
log_lh_time = toc(log_lh_tic);

% Sum up local cluster log-likelihoods
if (F.gossip)
    [log_lh_cluster, aggregate_error_ratio] = computeAggregateGossip(log_lh_cluster_ss, F.A, F.max_gossip_iter);
else
    log_lh_cluster = sum(log_lh_cluster_ss, 1);
end

% Now compute individual particle joint log-likelihood
% Construct the K-nearest graph for all the particles
% We actually find the k+1 nearest neighbors since particle i is the
% closest neighbor to particle i itself with 0 distance
% We thus ignore the first column of idx 
graph_tic = tic;
A = DelaunayGraph(x_predicted(1:2,:)');
% Construct Laplacian matrix
L = diag(sum(A,2)) - A;
graph_time = toc(graph_tic);

gamma_tic = tic;
% Solve the convex optimization
options =  optimoptions(@quadprog, 'Display','off');
gamma_approx = quadprog(L,[],[],[],C,log_lh_cluster',[], [], [], options)';
gamma_time = toc(gamma_tic);

gamma_exact = sum(log_lh_ss);
gamma_dif = gamma_approx-gamma_exact;

gamma = gamma_approx-max(gamma_approx);

try
    % Compute unnormalized posterior weight
    particle_weights = exp(gamma).*x_predicted(d+1,:);
catch ME
    save('error.mat');
end

% Give all particles equal weights if all particles have zero weight
if (sum(particle_weights) == 0)
    % This should never happen
    warning('All particle weights vanished');
    particle_weights = ones(1,N);
end

% Normalize the weights
particle_weights = particle_weights./sum(particle_weights); 

gamma_exact = gamma_exact - max(gamma_exact);
weight_exact = exp(gamma_exact).*x_predicted(d+1,:);
weight_exact = weight_exact/sum(weight_exact);

weight_dif = weight_exact-particle_weights;