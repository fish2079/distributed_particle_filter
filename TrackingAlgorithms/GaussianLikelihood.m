function [particle_weights, gamma, gamma_ss] = GaussianLikelihood(x_predicted, F, D, obs)
%   Function to compute the posterior weights of particles
%   The weight is computed exactly in a centralized manner
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
% Nov. 9th, 2017

d = size(x_predicted,1)-1;
gamma = zeros(1, size(x_predicted,2));

% Loop through each measurement
for j = 1:size(D.measurements,2)
    z_received = D.measurements(:,j);
    % Compute expected measurement
    z_expected = obs.model(x_predicted(1:d,:), D.sensorLoc(:,j), obs);
    
    % Compute the Gaussian log-likelihood
    z_dif = F.minus(z_received, z_expected);
    
    gamma = gamma + log(mvnpdf(z_dif', obs.mu', obs.R))';
    gamma_ss(j,:) = log(mvnpdf(z_dif', obs.mu', obs.R))';
end

gamma = gamma - max(gamma);

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