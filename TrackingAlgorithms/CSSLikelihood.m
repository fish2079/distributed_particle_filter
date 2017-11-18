function particle_weights = CSSLikelihood(x_predicted, F, D, obs)
%   Function to compute the approximate posterior particles weights
%   The log-likelihood is computed in a distributed manner using CSS
%   methods
%   Note that this function is designed for bearing-only tracking
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

% Pre-allocate a vector to store the new variance
variance_theta_dist = zeros(1, size(z,2));
R = zeros(1, size(z,2)); 

% Loop through each measurement
for i = 1:size(z,2)   
    variance_theta_dist(i) = sigma_theta^2;
    
    % distance between sensor and particles
    particle_dist = bsxfun(@minus, x_predicted(1:2,:), x_o(:,i));
    particle_dist_squared = particle_dist.^2;
%     E = sum(sum(particle_dist_squared))/size(x_predicted,2);
    E =  sum(particle_dist_squared,1)*x_predicted(5,:)';
    % new measurement variance
    R(i) = E*0.5*(1-exp(-2*variance_theta_dist(i))); 
end

% Compute local CSS for each sensor
% CSS_matrix is nb_sensors x 6 where each row contains the local CSS of
% one sensor
CSS_matrix = calculate_CSS(z(1,:), x_o, R);

% Consensus to compute the sum of local_log_likelihoods (multiplied by 
% number of nodes K to compute sum instead of average)
%[css_gossip, scalars1(1)] = gossip((CSS_matrix_full')*K, F);
%[css_consensus, scalars1(2)] = max_consensus_regular(css_gossip, F);

% Since we assume loseless communication, we can simply run a large number
% of gossip iterations and consensus iterations so that all sensors can
% obtain the correct sum
% Therefore, we skip the distributed computation and compute the sum
% exactly
CSS = sum(CSS_matrix,1)';

% once consensus is done, every sensor just calculates the particle weights
multi_matrix = [ones([size(x_predicted,2),1]), x_predicted(1,:)'.^2, x_predicted(2,:)'.^2,-2*x_predicted(1,:)'.*x_predicted(2,:)',2*x_predicted(1,:)',-2*x_predicted(2,:)'];
temp = -0.5*multi_matrix*CSS;

gamma = temp';

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

function CSS  = calculate_CSS(y, x_o, variance_theta)
% function to calculate the constraint sufficient statitsics for the
% distributed particle using the measurement and measurement noise
x_o = x_o(1:2,:);
Z = bsxfun(@times,[-cos(y);sin(y)],x_o);
Z = sum(Z,1);

G1 = Z.^2./variance_theta;
G2 = (cos(y)).^2./variance_theta;
G3 = (sin(y)).^2./variance_theta;
G4 = sin(y).*cos(y)./variance_theta;
G5 = Z.*cos(y)./variance_theta;
G6 = Z.*sin(y)./variance_theta;

CSS = [G1; G2; G3; G4; G5; G6;]';

