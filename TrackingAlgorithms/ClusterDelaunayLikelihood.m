function [particle_weights, gamma_dif, weight_dif, cluster_time, log_lh_time, graph_time, gamma_time, aggregate_error_ratio, errorNorm, smoothness, mth_eigvalue] = ClusterDelaunayLikelihood(x_predicted, F, D, obs)
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
    
    log_lh_ss(i,:) = log(mvnpdf(z_dif', obs.mu', obs.R)+realmin)';
  
    log_lh_cluster_ss(i,:) = C*log_lh_ss(i,:)';
end
log_lh_time = toc(log_lh_tic);

% Sum up local cluster log-likelihoods
if (F.gossip)
    [log_lh_cluster, aggregate_error_ratio] = computeAggregateGossip(log_lh_cluster_ss, F.A, F.max_gossip_iter);
else
    log_lh_cluster = sum(log_lh_cluster_ss, 1);
    % Inject perturbation in the results
    log_lh_cluster = log_lh_cluster + log_lh_cluster.*((rand(1,numel(log_lh_cluster))<0.5)*2-1)*F.perturbation;
    aggregate_error_ratio = zeros(1, F.cluster.k);
end

% Construct the KNN or Delaunchy triangulation graph for all the particles
% Construct the KNN or Delaunchy triangulation graph for all the particles
graph_tic = tic;
switch F.LA.graphMethod
    case 'KNN'
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
    case 'Delaunay'
        A = DelaunayGraph(x_predicted(1:2,:)');
    case 'Epsilon'
        A = EpsilonGraph(x_predicted(1:2,:)', F.LA.epsilon);
end

smoothness = sqrt(sum(log_lh_ss,1)*L*(sum(log_lh_ss,1)'));
mth_eigvalue = temp(F.LA.m,F.LA.m);

% Change to weighted adjacency matrix if needed
if (F.cluster.weightedEdge)
    for i=1:F.N
        x1 = x_predicted(1:2,i);
        x2 = x_predicted(1:2,A(i,:)>0);
        dist = sqrt((x1(1,:)-x2(1,:)).^2+(x1(2,:)-x2(2,:)).^2);
        switch F.cluster.weightedEdgeStyle
            case 1 % 1/dist
                A(i,A(i,:)>0) = 1./dist;
            case 2 % 1/dist^2
                A(i,A(i,:)>0) = 1./(dist.^2);
            case 3 % exp(-dist^2/sigma)
                A(i,A(i,:)>0) = exp(-dist.^2);
        end

    end
end

% Construct Laplacian matrix
L = diag(sum(A,2)) - A;
graph_time = toc(graph_tic);

gamma_tic = tic;
% Solve the convex optimization
options =  optimoptions(@quadprog, 'Display','off');
gamma_approx = quadprog(L,[],[],[],C,log_lh_cluster',[], [], [], options)';
gamma_time = toc(gamma_tic);

gamma_exact = sum(log_lh_ss, 1);
gamma_dif = norm(gamma_approx-gamma_exact);

gamma_noGossipError = quadprog(L,[],[],[],C,sum(log_lh_cluster_ss, 1)',[], [], [], options)';

% errorNorm(1) = norm(C);
% errorNorm(2) = norm(log_lh_cluster-sum(log_lh_cluster_ss, 1));
% yo = gamma_approx-gamma_noGossipError;
% errorNorm(6) = norm(yo);
% yo = yo-max(yo);
% errorNorm(7) = norm(yo);
% yo = exp(yo);
% errorNorm(8) = norm(yo);
% yo = 1-yo;
% errorNorm(3) = norm(yo);
% errorNorm(4) = mean(abs(yo));
% Psi = C';
% errorNorm(5) = norm(Psi/(Psi'*Psi));

llh_matrix_un = [gamma_noGossipError', gamma_approx', sum(log_lh_ss,1)'];
llh_matrix = bsxfun(@minus, llh_matrix_un, max(llh_matrix_un,[],1));
lh_matrix = exp(llh_matrix);
lh_matrix = bsxfun(@rdivide, lh_matrix, sum(lh_matrix,1));
llh_matrix = log(lh_matrix+realmin);
C_norm = llh_matrix_un - llh_matrix;

delta_m =  bsxfun(@rdivide,(llh_matrix(:,1)-llh_matrix(:,3)), llh_matrix(:,3));
delta_m(isinf(abs(delta_m)))=0;
errorNorm(1) = max(abs(delta_m));
delta_gossip = bsxfun(@rdivide, (llh_matrix(:,2)-llh_matrix(:,1)),llh_matrix(:,1));
delta_gossip(isinf(abs(delta_gossip)))=0;
errorNorm(2) = max(abs(delta_gossip));

errorNorm(3) = 0; 
errorNorm(4) = 0; 
delta_llh = (llh_matrix(:,2)-llh_matrix(:,3))./(llh_matrix(:,3));
delta_llh(isinf(delta_llh))=0;
errorNorm(5) = max(abs(delta_llh));
errorNorm(6) = (1+errorNorm(1))*errorNorm(2)+errorNorm(1);

errorNorm(7) = min(sqrt(exp(llh_matrix(:,1)-llh_matrix(:,2))));
errorNorm(8:13) = 0;

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

gamma_exact = gamma_exact - max(gamma_exact);
weight_exact = exp(gamma_exact).*x_predicted(d+1,:);
weight_exact = weight_exact/sum(weight_exact);

weight_dif = norm(weight_exact-particle_weights);