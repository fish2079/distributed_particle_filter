function [particle_weights, aggregate_error_ratio, errorNorm] = CSSLikelihood(x_predicted, F, D, obs)
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

for i=1:numel(D.sensorID)    
    z_received = D.measurements(:,i);
    % Compute expected measurement
    z_expected = obs.model(x_predicted(1:d,:), D.sensorLoc(:,i), obs);
    
    % Compute the Gaussian log-likelihood
    z_dif = F.minus(z_received, z_expected);
    
    log_lh_ss(i,:) = log(mvnpdf(z_dif', obs.mu', obs.R)+realmin)';
end

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
if (F.gossip)
    [CSS, aggregate_error_ratio] = computeAggregateGossip(CSS_matrix, F.A, F.max_gossip_iter);
    CSS = CSS';
else
    CSS = sum(CSS_matrix,1)';
    % Inject perturbation in the results
    CSS = CSS + CSS.*((rand(6,1)<0.5)*2-1)*F.perturbation;
    aggregate_error_ratio = zeros(1,6);
end

% once consensus is done, every sensor just calculates the particle weights
multi_matrix = [ones([size(x_predicted,2),1]), x_predicted(1,:)'.^2, x_predicted(2,:)'.^2,-2*x_predicted(1,:)'.*x_predicted(2,:)',2*x_predicted(1,:)',-2*x_predicted(2,:)'];
temp = -0.5*multi_matrix*CSS;

% errorNorm(1) = norm(-0.5*multi_matrix);
% errorNorm(2) = norm(CSS-sum(CSS_matrix,1)');
% yo = -0.5*multi_matrix*(CSS-sum(CSS_matrix,1)');
% errorNorm(6) = norm(yo);
% yo = yo-max(yo);
% errorNorm(7) = norm(yo);
% yo = exp(yo);
% errorNorm(8) = norm(yo);
% yo = 1-yo;
% errorNorm(3) = norm(yo);
% errorNorm(4) = mean(abs(yo));
% Psi = -0.5*multi_matrix;
% errorNorm(5) = norm(Psi/(Psi'*Psi));

Psi = -0.5*multi_matrix;
alpha_true = sum(CSS_matrix,1)';
alpha_gossip = CSS;
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

