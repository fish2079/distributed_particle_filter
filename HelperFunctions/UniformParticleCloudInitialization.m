function x_new = UniformParticleCloudInitialization(F)
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

x_new = [F.initial_upper_limit-F.initial_lower_limit].*rand(F.N,F.d)+F.initial_lower_limit;
x_new = x_new';
x_new(F.d+1,:) = 1/F.N;
