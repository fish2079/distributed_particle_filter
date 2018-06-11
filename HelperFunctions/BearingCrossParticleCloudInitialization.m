function x_new = BearingCrossParticleCloudInitialization(F)
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

delta = [sin(F.initial_meas)*200;cos(F.initial_meas)*200]; 

sensor_max = F.initial_sensor+delta;

nb_sensors = size(sensor_max,2);

x_intersect = [];
for i=nb_sensors:-1:1
    for j=i-1:-1:1
        x_intersect = [x_intersect, InterX([F.initial_sensor(:,i),sensor_max(:,i)],[F.initial_sensor(:,j),sensor_max(:,j)])];
    end
end

[mu(1),sigma(1)]=normfit(x_intersect(1,:));
[mu(2),sigma(2)]=normfit(x_intersect(2,:));

x_new = mvnrnd([mu(1),mu(2), 0,0]',diag([sigma(1), sigma(2), 5,5].^2),F.N)';
x_new = x_new + F.sigma_noise*randn(F.d,F.N);


%x_new = mvnrnd(F.initial_mu, F.initial_R, F.N)';
x_new(F.d+1,:) = 1/F.N;
