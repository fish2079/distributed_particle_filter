function x_new = GaussianParticleCloudInitialization(F)
%   Function to initialize the initial particle cloud based on
%   provided algorithm parameters
%
%   Inputs:
%       F: struct, contain all parameters needed for initialization    
%
%   Outputs:
%       x_new: (d+1)xN matrix of initial particles, columns 1 to d are the 
%       particles and column 5 is the particle weight
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

x_new = mvnrnd(F.initial_mu, F.initial_R, F.N)';
x_new(F.d+1,:) = 1/F.N;
