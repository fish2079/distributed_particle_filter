function [x_new] = propagate_cv_ct(x_old, dynamic)
%   Function to implement the dynamic model of maneuvering target at
%   constant velocity with ppossibility of clockwise turn
%
%   Input:
%       x_old: 4xN matrix of current target states where N is the number of
%       particles
%       dynamic: struct containing all relevant parameters of target dynamic
%       model
%
%   Output: 
%       x_new: 4xN matrix of new target states after propagation
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

% Load the model parameters
a = dynamic.a; % turning rate
sigma_a = dynamic.sigma_a; % system noise
p = dynamic.p; % probability of turning
T = dynamic.T; % sampling interval, minutes

d = size(x_old,1); % target state dimension
N = size(x_old,2); % Number of particles

% Initialize output
x_new = zeros(d, N);

% Generate covariance matrix and process noise
covariance = [T^3/3, T^2/2;T^2/2, T]*sigma_a^2;
ww1 = mvnrnd([0,0,], covariance, N)'; 
ww2 = mvnrnd([0,0,], covariance, N)'; 
ww = [ww1(1,:);ww2(1,:);ww1(2,:);ww2(2,:)];

% Flip coins to see which particles are turning
turning = (rand(N,1) <= p);
notturning = ~turning;

% Update the non-turning particles
% Compute the dispacements  for non-turning particles
displacement(1:2, notturning)=[T.*x_old(3, notturning);T.*x_old(4, notturning)];
% For non-turning particles, their velocity does not change
x_new(3, notturning) = x_old(3, notturning);
x_new(4, notturning) = x_old(4, notturning);

% Update the turning particles
omega = a./sqrt(sum(x_old(3:4,turning).^2,1)); % turning rate
sinOmegaT = sin(omega.*T);
cosOmegaT = cos(omega.*T);
tempa = ones(size(omega));
tempb = ones(size(omega));
tempa(omega ~= 0) = sinOmegaT(omega ~= 0)./omega(omega ~= 0);
tempb(omega ~= 0) = (1 - cosOmegaT(omega ~= 0))./omega(omega ~= 0);
% Compute the dipslacements for turning particles
displacement(1:2, turning)=[tempa.*x_old(3,turning) - tempb.*x_old(4,turning);tempb.*x_old(3,turning) + tempa.*x_old(4,turning)];
% Compute the new velocity for turning particles
x_new(3,turning) = cosOmegaT.*x_old(3,turning) - sinOmegaT.*x_old(4,turning);
x_new(4,turning) = sinOmegaT.*x_old(3,turning) + cosOmegaT.*x_old(4,turning);
% Compute the new position for turning particles
x_new(1,:)=x_old(1,:)+displacement(1,:);
x_new(2,:)=x_old(2,:)+displacement(2,:);

% Inject process noise
x_new = x_new + ww;

