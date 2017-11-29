function [particle_weights, cluster_time, KNN_time, gamma_time] = ClusterLikelihood(x_predicted, F, D, obs)
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
for i=1:numel(D.sensorID)
    D_single.measurements = D.measurements(:,i);
    D_single.sensorID = D.sensorID(i);
    D_single.sensorLoc = D.sensorLoc(:,i);
    log_lh_ss(i,:) = log(GaussianLikelihood(x_predicted, F, D_single, obs)+realmin);
    log_lh_cluster_ss(i,:) = C*log_lh_ss(i,:)';
end

% Sum up local cluster log-likelihoods
log_lh_cluster = sum(log_lh_cluster_ss, 1);

% Now compute individual particle joint log-likelihood
% Construct the K-nearest graph for all the particles
% We actually find the k+1 nearest neighbors since particle i is the
% closest neighbor to particle i itself with 0 distance
% We thus ignore the first column of idx 
KNN_tic = tic;
idx = knnsearch(x_predicted(1:2,:)', x_predicted(1:2,:)','k',F.cluster.KNN+1);
idx = idx(:,2:end);
KNN_time = toc(KNN_tic);

% Now construct the adjacency matrix
A = zeros(F.N, F.N);
for i=1:F.N
    % particle i is connected to its K nearest neighbor
    A(i,idx(i,:)) = 1;
    % the connection is symmetric
    A(idx(i,:),i) = 1;
end

% Construct Laplacian matrix
L = diag(sum(A,2)) - A;

gamma_tic = tic;
% Solve the convex optimization
options =  optimoptions(@quadprog, 'Display','off');
gamma = quadprog(L,[],[],[],C,log_lh_cluster',[], [], [], options)';
gamma_time = toc(gamma_tic);

gamma = gamma-max(gamma);

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