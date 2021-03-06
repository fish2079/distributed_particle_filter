function particle_weights = CSSLikelihoodExact(x_predicted, F, D, obs)
%   Function to compute the approximate posterior particles weights
%   The log-likelihood is computed in a centralized manner using CSS
%   methods
%   Note that this function is designed for bearing-only tracking
%   In addition, we compute the per-particle variance rather than using
%   weighted average
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
% Nov. 9th, 2017

d = size(x_predicted,1)-1;

z = D.measurements;

% Extra parameters
sigma_theta = sqrt(obs.R);
% sigma_a = dynamic.sigma_a;
x_o = D.sensorLoc;

gamma = zeros(1, size(x_predicted,2));
for i= 1:size(z,2)
    % Compute transformed received measurement
    Z(i) = sum(bsxfun(@times,[-cos(z(i));sin(z(i))], x_o(:,i)),1);

    % Compute transformed expected measurement
    Z_exp(i,:) = sum(bsxfun(@times,[-cos(z(i));sin(z(i))],x_predicted(1:2,:)),1);
   
%     % distance between sensor and particles
    particle_dist = bsxfun(@minus, x_predicted(1:2,:), x_o(:,i));
%     particle_dist_squared = particle_dist.^2;

    % new per-particle measurement variance
    R(i,:) = sum(particle_dist.^2,1).*0.5*(1-exp(-2*sigma_theta^2)); 

    gamma = gamma-(bsxfun(@minus, Z(i), Z_exp(i,:))).^2/2./R(i,:)-log(sqrt(2*pi*R(i,:)));
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