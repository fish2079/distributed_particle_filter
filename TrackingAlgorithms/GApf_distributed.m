function [x_updated, details, errorNorm] = GApf_distributed(x_old, F, D, dynamic, obs, details)
%   This function implements one time step of the distributed Gaussian
%   approximation partcile filter
%
%   Inputs:
%       x_old: (d+1)*N matrix of particles from previous time step where N
%              is the number of particles and d is the dimension of the  
%              state. The d+1 row is the particle weights
%       F: Struct containing filter parameters
%       D: Struct containing measurement data
%       dynamics: Struct containing dynamic model parameters
%       obs: Struct containing measurement model paraleters
%
%   Outputs:
%       x_updated: (d+1)*N matrix of updated particles where N is the number of
%                  particles and d is the dimension of the state. The d+1 
%                  row is the updated particle weights
%
% Oreshkin and Coates, Fusion 2010. (as described in paper)

N = F.N; % number of particles
nb_sensors = size(obs.sensorPos,2); % number of sensors
d = F.d; % state dimension

% Generate predicted particles
x_predicted = dynamic.model(x_old(1:d,:), dynamic);

% Initial Gaussian distribution parameters for all sensors
MUini = zeros(d,nb_sensors);
Rini = zeros(d,d,nb_sensors);
Rini_inv = zeros(d,d,nb_sensors);

% Loop through each sensor
for i = 1:nb_sensors 
    % Compute expected measurement
    z_received = D.measurements(:,i);
    z_expected = obs.model(x_predicted(1:d,:), D.sensorLoc(:,i), obs);
    z_dif = F.minus(z_received, z_expected);

    % Compute the Gaussian log-likelihood
    gammai = log(mvnpdf(z_dif', obs.mu', obs.R)+realmin)';
    gammai = gammai - max(gammai); 
    pw = exp(gammai);
    
    pw = x_old(d+1,:) .* pw;
    pw = pw ./ sum(pw);
        
    % Sample auxiliary particles
    I = randsample((1:N)', N, true, pw);

    % Generate true predicted particles from auxiliary particles
    x_intermediate = dynamic.model(x_old(1:d, I), dynamic);

    % Step 2, compute weights for predicted particles
    z_expected = obs.model(x_intermediate(1:d,:), D.sensorLoc(:,i), obs);
    z_dif = F.minus(z_received, z_expected);
    log_lh_ss(i,:) = log(mvnpdf(z_dif', obs.mu', obs.R)+realmin)';
    log_lh_ss(i,:) = bsxfun(@minus, log_lh_ss(i,:), max(log_lh_ss(i,:)));
    
    lh_ss(i,:) = exp(log_lh_ss(i,:));
    lh_ss(i,:) = bsxfun(@rdivide, lh_ss(i,:), pw(I));
    lh_ss(i,:) = lh_ss(i,:)/sum(lh_ss(i,:));
    
    % Compute Gaussian approximaiton parameters
    [MUini(:,i), Rini(:,:,i)] = save_gs_stat([x_intermediate; lh_ss(i,:)]');
    Rini_inv(:,:,i) = chol_inv(Rini(:,:,i));
end

% Use gossiping and max consensus to fuse the Gaussian parameters
for i = 1:nb_sensors
    Rave_matrix(:,:,i) = Rini_inv(:,:,i);
    MUave_matrix(:,i) = Rini_inv(:,:,i) * MUini(:,i);
    coefficients_matrix(:,i) = [MUave_matrix(:,i);reshape(Rave_matrix(:,:,i),d^2,1)];
end

if (F.gossip)
    [coefficients, aggregate_error_ratio] = computeAggregateGossip(coefficients_matrix', F.A, F.max_gossip_iter);
    if (isfield(details,'aggregate_error_ratio'))
        details.aggregate_error_ratio = [details.aggregate_error_ratio; aggregate_error_ratio];
    else
        details.aggregate_error_ratio = aggregate_error_ratio;
    end
else
    coefficients = sum(coefficients_matrix,2)';
    aggregate_error_ratio = zeros(1, size(coefficients,1));
end

% Construct the aggreate Gaussian distribution parameters
MU_ave = coefficients(:,1:d)';
R_ave = reshape(coefficients(:,d+1:end),d,d)';
R_ave = chol_inv(R_ave);
MU_ave = R_ave*MU_ave;
R_ave = nb_sensors*R_ave;
   
% Perturbation noise to inject into the particles
noise1 = F.sigma_noise*randn(N,d);
x_updated = [mvnrnd(MU_ave, R_ave, N)+noise1, ones(N,1)/N]';

% Debug codes
coefficients_exact = sum(coefficients_matrix,2)';
MU_ave_exact = coefficients_exact(:,1:d)';
R_ave_exact = reshape(coefficients_exact(:,d+1:end),d,d)';
R_ave_exact = chol_inv(R_ave_exact);
MU_ave_exact = R_ave_exact*MU_ave_exact;
R_ave_exact = nb_sensors*R_ave_exact;

gamma_noGossipError = mvnpdf(x_intermediate', MU_ave_exact',R_ave_exact);
gamma_approx = mvnpdf(x_intermediate', MU_ave',R_ave);

llh_matrix_un = [gamma_noGossipError, gamma_approx, sum(log_lh_ss,1)'];
llh_matrix = bsxfun(@minus, llh_matrix_un, max(llh_matrix_un,[],1));
lh_matrix = exp(llh_matrix);
lh_matrix = bsxfun(@rdivide, lh_matrix, sum(lh_matrix,1));
llh_matrix = log(lh_matrix+realmin);
% C_norm = llh_matrix_un - llh_matrix;

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
errorNorm(8:12) = 0;

gamma_dif = norm(llh_matrix(:,3)-llh_matrix(:,2));
if (isfield(details,'gamma_dif'))
    details.gamma_dif = [details.gamma_dif, gamma_dif];
else
    details.gamma_dif = gamma_dif;
end

weight_dif = norm(lh_matrix(:,3)-lh_matrix(:,2));
if (isfield(details,'weight_dif'))
    details.weight_dif = [details.weight_dif, weight_dif];
else
    details.weight_dif = weight_dif;
end