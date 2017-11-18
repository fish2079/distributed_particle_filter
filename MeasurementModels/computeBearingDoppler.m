function measurements = computeBearingDoppler(x_t, sensorPos, obs)
%   Measurement function to compute the bearing and Doppler frequency shift 
%   from target(s) to sensor(s)
%
%   Inputs:
%       x_t: 4-by-N matrix of particle/target states
%       sensorPos: 2-by-nb_sensors matrix of sensor positions
%       obs: struct containing all measurement model parameters
%
%   Outputs:
%       measurements: 2nb_sensors x nb_particles matrix of measurements
%       Each odd row contains bearing measurements (in radians). Each even 
%       row contains Doppler measurements. Every 2 row corresponds to one 
%       single sensor' measurements. 
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

N = size(x_t,2); % Number of particles/targets
S = size(sensorPos,2); % Number of sensors
measurements = [];

% Loop through each sensors
for kk=1:S
    % Compute sensor kk's bearings relative to all targets
    difX = x_t(1,:) - sensorPos(1,kk);
    difY = x_t(2,:) - sensorPos(2,kk);
    range = sqrt(difX.^2+difY.^2);
    bearing = atan2(difX, difY);
    doppler = 2./obs.lambda.*((difX.*x_t(3,:)+difY.*x_t(4,:))./range);
    % If the target and the sensor have 0 distance, then replace NaN by 0
    % doppler shift
    doppler(range==0) = 0;
    measurements = [measurements; bearing; doppler];
end