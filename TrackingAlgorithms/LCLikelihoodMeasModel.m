function [particle_weights, Hx_ss_dif, aggregate_error_ratio] = LCLikelihoodMeasModel(x_predicted, F, D, obs)
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

z = D.measurements;
d = size(x_predicted,1)-1;

% Precompute inverse of measurement noise covariance matrix
R_inv = inv(obs.R);

% Construct the Psi matrix
degree_matrix = combinator(F.LC.max_degree+1,2,'p','r')'-1;
for i=1:size(degree_matrix,2)
    Psi(:,i) = prod(bsxfun(@power, x_predicted(1:2,:), degree_matrix(:,i)),1);
end
% Psi = exp(log(x_predicted(1:2,:))'*degree_matrix);
% for i=1:size(F.LC.basis,2)
%     
%     Psi(:,i) = mvnpdf(x_predicted(1:2,:)', F.LC.basis(:,i)', F.LC.R);
% end

% Compute the coefficients for likelihood consensus
for i=1:size(z,2)
    % Compute the true value beta vector
    % Since we are approximating the measurement model, beta contains the
    % expected measurements
    beta_LC(:,:,i) = obs.model(x_predicted(1:d,:), D.sensorLoc(:,i), obs);
    
    % Compute the coefficients
    alpha_LC(:,:,i) = mldivide(Psi,beta_LC(:,:,i)');
end   

% Compute the first local statistics z'Q^{-1}\alpha
for i=1:size(z,2)
    % Compute z'Q^{-1}\alpha
    zHx_ss(:,i) = z(:,i)'*R_inv*alpha_LC(:,:,i)';
    Hx_ss_est(:,:,i) = (Psi*alpha_LC(:,:,i))';
    % This line of debug code is specific to bearing model
    Hx_ss_dif(:,i) = wrapToPi(beta_LC(:,:,i)-Hx_ss_est(:,:,i));
end

% Compute the second local statistics \alpha'Q^{-1}\alpha
% and corresponding modified basis funciton
% idx = 1;
% for i=1:size(F.LC.basis,2)
%     for j=1:size(F.LC.basis,2)
%         for zz=1:size(z,2)
%             HxHx_ss(zz,idx) = -0.5*alpha_LC(i,:,zz)*R_inv*alpha_LC(j,:,zz)';          
%         end 
%         Psi_HxHx(:,idx) = Psi(:,i).*Psi(:,j);
%         idx = idx+1;
%     end
% end

% Compute the second local statistics \alpha'Q^{-1}\alpha
HxHx_ss = [];
[a,b]=meshgrid(1:size(Psi,2),1:size(Psi,2));

Psi_extended = Psi(:,a(:)).*Psi(:,b(:));
for zz=1:size(z,2)
    alpha_LC_extLeft = alpha_LC(a(:),:,zz);
    alpha_LC_extRight = alpha_LC(b(:),:,zz);
    temp = -0.5*alpha_LC_extLeft*R_inv;
    temp = temp'.*alpha_LC_extRight';
    HxHx_ss(zz,:) =  sum(temp,1);
end

if (F.gossip)
    [zHx, aggregate_error_ratio1] = computeAggregateGossip(zHx_ss', F.A, F.max_gossip_iter);
    [HxHx, aggregate_error_ratio2] = computeAggregateGossip(HxHx_ss, F.A, F.max_gossip_iter);
    zHx = zHx';
    HxHx = HxHx';
    aggregate_error_ratio = [aggregate_error_ratio1, aggregate_error_ratio2];
else
    zHx = sum(zHx_ss,2);
    HxHx = sum(HxHx_ss,1)';
     % Inject perturbation in the results
    zHx = zHx + zHx.*((rand(size(zHx,1),1)<0.5)*2-1)*F.perturbation;
    HxHx = HxHx + HxHx.*((rand(size(HxHx,1),1)<0.5)*2-1)*F.perturbation;
    
    aggregate_error_ratio = zeros(1, numel(zHx)+numel(HxHx));
end

gamma = Psi*zHx + Psi_extended*HxHx;

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
