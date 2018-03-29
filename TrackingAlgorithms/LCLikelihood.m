function [particle_weights, aggregate_error_ratio] = LCLikelihood(x_predicted, F, D, obs)
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
