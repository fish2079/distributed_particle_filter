function range = computeRange(x_t, sensorPos, obs)
%   Function to compute the range from target(s) to sensor(s)
%
%   Inputs:
%       x_t: 4-by-N matrix of particle/target states
%       sensorPos: 2-by-nb_sensors matrix of sensor positions
%       obs: struct containing all measurement model parameters
%
%   Outputs:
%       bearing: nb_sensors x N matrix of bearing measurements (in radians)
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

N = size(x_t,2); % Number of particles/targets
S = size(sensorPos,2); % Number of sensors
range = zeros(S,N);

% Loop through each sensors
for kk=1:S
    % Compute sensor kk's bearings relative to all targets
    difX = x_t(1,:) - sensorPos(1,kk);
    difY = x_t(2,:) - sensorPos(2,kk);
    range(kk,:) = sqrt(difX.^2+difY.^2);
end


