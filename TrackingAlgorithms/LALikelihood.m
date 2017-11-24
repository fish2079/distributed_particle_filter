function particle_weights = LALikelihood(x_predicted, F, D, obs)
%   Function to compute the approximate posterior particles weights
%   The log-likelihood is computed in a distributed manner using Laplacian
%   approximation methods
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

% Start by constructing the K-nearest graph for all the particles
% We actually find the k+1 nearest neighbors since particle i is the
% closest neighbor to particle i itself with 0 distance
% We thus ignore the first column of idx 
idx = knnsearch(x_predicted(1:2,:)', x_predicted(1:2,:)','k',F.LA.KNN+1);
idx = idx(:,2:end);

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

% Do eigenvalue decomposition of Laplacian matrix
[V_full,eigValues] = eig(L);

% Select the m smallest eigenvectors;
V = V_full(:,1:F.LA.m);

% Now have each sensor compute local log-likelihood using only local
% measurements
for i=1:numel(D.sensorID)
    D_single.measurements = D.measurements(:,i);
    D_single.sensorID = D.sensorID(i);
    D_single.sensorLoc = D.sensorLoc(:,i);
    log_lh_ss(i,:) = log(GaussianLikelihood(x_predicted, F, D_single, obs));
end

% Compute local coefficients
alpha_ss = V'*log_lh_ss';

% Sum up the local coefficients
alpha = sum(alpha_ss,2);

% Compute approximate global joint log-likelihood
gamma = (V*alpha)';

gamma = gamma-max(gamma);

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