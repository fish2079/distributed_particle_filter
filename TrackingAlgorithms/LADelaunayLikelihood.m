function [particle_weights, gamma_dif, weight_dif, log_lh_time, graph_time, eig_time, aggregate_error_ratio, errorNorm] = LADelaunayLikelihood(x_predicted, F, D, obs)
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

% First have each sensor compute local log-likelihood using only local
% measurements
log_lh_tic = tic;
for i=1:numel(D.sensorID)
%     D_single.measurements = D.measurements(:,i);
%     D_single.sensorID = D.sensorID(i);
%     D_single.sensorLoc = D.sensorLoc(:,i);
%     log_lh_ss_approx(i,:) = log(GaussianLikelihood(x_predicted, F, D_single, obs)+realmin);
    
    z_received = D.measurements(:,i);
    % Compute expected measurement
    z_expected = obs.model(x_predicted(1:d,:), D.sensorLoc(:,i), obs);
    
    % Compute the Gaussian log-likelihood
    z_dif = F.minus(z_received, z_expected);
    
    log_lh_ss(i,:) = log(mvnpdf(z_dif', obs.mu', obs.R)+realmin)';
end

log_lh_time = toc(log_lh_tic);

% Construct the KNN or Delaunchy triangulation graph for all the particles
graph_tic = tic;
if (F.LA.KNNgraph)
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
else
    A = DelaunayGraph(x_predicted(1:2,:)');
end

% Change to weighted adjacency matrix if needed
if (F.LA.weightedEdge)
    for i=1:F.N
        x1 = x_predicted(1:2,i);
        x2 = x_predicted(1:2,A(i,:)>0);
        dist = sqrt((x1(1,:)-x2(1,:)).^2+(x1(2,:)-x2(2,:)).^2);
        A(i,A(i,:)>0) = 1./dist;
    end
end

% Construct Laplacian matrix
L = diag(sum(A,2)) - A;
graph_time = toc(graph_tic);

% Do eigenvalue decomposition of Laplacian matrix
eig_time_tic = tic;
[V_full,~] = eig(L);
% V_full = F.LA.V_full;
eig_time = toc(eig_time_tic);

% Select the m smallest eigenvectors;
V = V_full(:,1:F.LA.m);

% Compute local coefficients
alpha_ss = V'*log_lh_ss';

% Sum up the local coefficients
if (F.gossip)
    [alpha, aggregate_error_ratio] = computeAggregateGossip(alpha_ss', F.A, F.max_gossip_iter);
    alpha = alpha';
else
    alpha = sum(alpha_ss,2);
    % Inject perturbation in the results
    alpha = alpha + alpha.*((rand(numel(alpha),1)<0.5)*2-1)*F.perturbation;
    aggregate_error_ratio = zeros(1,F.LA.m);
end

Psi = V;
alpha_true = sum(alpha_ss,2);
alpha_gossip = alpha;
llh_matrix = [Psi*alpha_true, Psi*alpha_gossip, sum(log_lh_ss,1)'];
lh_matrix = exp(llh_matrix);
lh_matrix = lh_matrix./sum(lh_matrix,1);
llh_matrix = log(lh_matrix+realmin);
% delta_m=(Psi*alpha_true-sum(log_lh_ss,1)')./sum(log_lh_ss,1)';
delta_m = (llh_matrix(:,1)-llh_matrix(:,3))./llh_matrix(:,3);
delta_m(isinf(abs(delta_m)))=0;
% delta_gossip=(Psi*alpha_gossip-Psi*alpha_true)./(Psi*alpha_true);
delta_gossip = (llh_matrix(:,2)-llh_matrix(:,1))./llh_matrix(:,1);
delta_gossip(isinf(abs(delta_gossip)))=0;
errorNorm(1) = max(abs(delta_m));
errorNorm(2) = max(abs(delta_gossip));
for i=1:numel(delta_m)
%     tempUpper(i) = norm(Psi(i,:)')*norm(Psi(i,:))*norm(alpha_gossip-alpha_true)/norm(Psi(i,:)'*Psi(i,:)*alpha_true);
    tempUpper(i) = norm(Psi(i,:)')*norm(Psi(i,:))*norm(alpha_gossip-alpha_true)/norm(Psi(i,:)'*llh_matrix(i,1));
%     tempLower(i) = norm(Psi(i,:)'*Psi(i,:)*alpha_gossip-Psi(i,:)'*Psi(i,:)*alpha_true)/norm(Psi(i,:)')/norm(Psi(i,:)*alpha_true);
    tempLower(i) = norm(Psi(i,:)'*llh_matrix(i,1)-Psi(i,:)'*llh_matrix(i,2))/norm(Psi(i,:)')/norm(llh_matrix(i,1));
end
tempUpper(isinf(abs(tempUpper)))=0;
tempLower(isinf(abs(tempLower)))=0;
errorNorm(3) = max(abs(tempUpper));
errorNorm(4) = max(abs(tempLower));
delta_llh = (llh_matrix(:,2)-llh_matrix(:,3))./(llh_matrix(:,3));
delta_llh(isinf(delta_llh))=0;
errorNorm(5) = max(abs(delta_llh));
errorNorm(6) = (1+errorNorm(1))*errorNorm(2)+errorNorm(1);

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

gamma_dif = norm(gamma_approx-gamma_exact);

gamma_exact = gamma_exact-max(gamma_exact);
weight_exact = exp(gamma_exact).*x_predicted(d+1,:);
weight_exact = weight_exact/sum(weight_exact);

weight_dif = norm(weight_exact-particle_weights);