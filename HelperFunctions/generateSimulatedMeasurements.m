function D = generateSimulatedMeasurements(S, obs)
%   Function to generate sensor measurements based on the given track and
%   measurement model/parameters
%   
%   Input:
%       S: struct containing the target track and sensor positions
%       obs: struct containing all measurement model parameters
%
%   Output:
%       D: struct containing the generated measurements and receiving
%       sensors
%
%   At each time step one measurement is generated per sensor. The
%   measurement depends on the specified model and parameters. 
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

d = size(obs.R,1);

% Allocate space for returned variables
D.measurements = cell(1,S.nb_steps);
D.sensorID = cell(1,S.nb_steps);

% Loop to generate data
for t=1:S.nb_steps
    % Number of measurements per sensor
    % For now, assume 1 measurement per sensor
    % We ignore missed detection and clutter measurement
    numMeas = ones(1, S.nb_sensors);
    
    % Allocate space
    D.measurements{t} = [];
    D.sensorID{t} = [];
    D.sensorLoc{t} = [];
    D.noise{t} = [];
    
    % Loop over sensors and generate measurements
    for s=1:S.nb_sensors
        % Measurements for sensor s
        if numMeas(s) > 0
            % compute noiseless measurement
            z_noiseless = obs.model(S.x_t(:,t), obs.sensorPos(:,s), obs);
            % compute Gaussian measurement noise
            noise = mvnrnd(obs.mu', obs.R, numMeas(s))';
            D.measurements{t} = [D.measurements{t}, bsxfun(@plus, z_noiseless, noise)];
            D.sensorID{t} = [D.sensorID{t}, s*ones(1,numMeas(s))];
            D.sensorLoc{t} = [D.sensorLoc{t}, repmat(obs.sensorPos(:,s), 1, numMeas(s))];
            D.noise{t} = [D.noise{t}, noise];
        end        
    end
end



