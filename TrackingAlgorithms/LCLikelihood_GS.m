function [particle_weights, aggregate_error_ratio, errorNorm] = LCLikelihood_GS(x_predicted, F, D, obs)
%   Function to compute the posterior weights of particles
%   The weight is computed using likelihood consensus to approximate the
%   measurement model
%
%   Inputs:
%       x_predicted: (d+1)-by-N matrix of particle states
%       F: Struct containing filter parameters
%       D: Struct containing measurement data
%       obs: Struct containing measurement model paraleters
%
% Output:
%       particle_weights: 1-by-N row vector of posterior particle weights
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 13th, 2017

% z = D.measurements;
d = size(x_predicted,1)-1;

% Precompute inverse of measurement noise covariance matrix
% R_inv = inv(obs.R);

% Construct the Psi matrix
degree_matrix = combinator(F.LC.max_degree+1,2,'p','r')'-1;
for i=1:size(degree_matrix,2)
    Psi(:,i) = prod(bsxfun(@power, x_predicted(1:2,:), degree_matrix(:,i)),1);
end
Psi = GramSchmidt(Psi);
% Psi_sum = sum(Psi,1);
% Psi_normalized = Psi./Psi_sum;

% Compute the local log-likelihoods
for i=1:numel(D.sensorID)  
    z_received = D.measurements(:,i);
    % Compute expected measurement
    z_expected = obs.model(x_predicted(1:d,:), D.sensorLoc(:,i), obs);
    
    % Compute the Gaussian log-likelihood
    z_dif = F.minus(z_received, z_expected);
    
    log_lh_ss(i,:) = log(mvnpdf(z_dif', obs.mu', obs.R)+realmin)';
    
%     alpha_LC(:,i) = mldivide(Psi_normalized,log_lh_ss(i,:)');
    alpha_LC(:,i) = mldivide(Psi,log_lh_ss(i,:)');
%     temp_alpha_LC(:,i) = mldivide((Psi'*Psi),Psi')*log_lh_ss(i,:)';
end

if (F.gossip)
    [alpha_LC_aggregate, aggregate_error_ratio] = computeAggregateGossip(alpha_LC', F.A, F.max_gossip_iter);
    alpha_LC_aggregate = alpha_LC_aggregate';
else
    alpha_LC_aggregate = sum(alpha_LC,2);
    % Inject perturbation in the results
    alpha_LC_aggregate = alpha_LC_aggregate + alpha_LC_aggregate.*((rand(size(alpha_LC_aggregate,1),1)<0.5)*2-1)*F.perturbation;
    aggregate_error_ratio = zeros(1, numel(alpha_LC_aggregate));
end

% errorNorm(1) = norm(Psi);
% errorNorm(2) = norm(alpha_LC_aggregate-sum(alpha_LC,2));
% yo = Psi*(alpha_LC_aggregate-sum(alpha_LC,2));
% errorNorm(6) = norm(yo);
% yo = yo-max(yo);
% errorNorm(7) = norm(yo);
% yo = exp(yo);
% errorNorm(8) = norm(yo);
% yo = 1-yo;
% errorNorm(3) = norm(yo);
% errorNorm(4) = mean(abs(yo));
% errorNorm(5) = norm(Psi/(Psi'*Psi));

alpha_true = sum(alpha_LC,2);
alpha_gossip = alpha_LC_aggregate;
llh_matrix = [Psi*alpha_true, Psi*alpha_gossip, sum(log_lh_ss,1)'];
llh_matrix = llh_matrix-max(llh_matrix);
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

gamma = Psi*alpha_LC_aggregate;
% gamma = Psi_normalized*alpha_LC_aggregate;
gamma = gamma' - max(gamma);

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
