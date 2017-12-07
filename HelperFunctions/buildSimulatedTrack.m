function S_output = buildSimulatedTrack(S, dynamic)
%   Function to generate a simulated target track based on the parameters
%   in the structs S, F, and dynamic
%   
%   Input:
%       S: struct containing all relevant parameters for the track
%       dynamic: struct containing all target transition model parameters
%
%   Output:
%       S_output: the same struct as input s with the additional field of
%       x_t
%
% Jun Ye Yu
% McGill University
% jun.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

S_output = S;
d = 4; % target state dimension

%% Construct tracking area
% The initial target position is at the center of tracking area
xMin = S.initial(1)-S.area_length*0.5;
xMax = S.initial(1)+S.area_length*0.5;

% yMin = S.initial(2)-S.area_length*0.5;
% yMax = S.initial(2)+S.area_length*0.5;

yMin = S.initial(2)-S.area_length*0.15;
yMax = S.initial(2)+S.area_length*0.85;

% yMin = S.initial(2)+S.area_length*1.75;
% yMax = S.initial(2)+S.area_length*2.75;

%% Place sensors in the tracking area
% Place sensors in a grid if needed
if (S.grid_sensor)
    rootK = sqrt(S.nb_sensors);
    sensorX = xMin:S.area_length/(rootK-1):xMax;
    sensorY = yMin:S.area_length/(rootK-1):yMax;   
    [sensorX, sensorY] = meshgrid(sensorX, sensorY);
    sensorPos = [sensorX(:), sensorY(:)]';
else
    % Alternatively, place sensors randomly in the tracking field
    offset = rand(S.nb_sensors, 2);
    sensorPos = [xMin + offset(:,1)*area_length, yMin + offset(:,2)*area_length]';
end
% Store all sensor positions
S_output.sensorPos = sensorPos;

%% Construct target track
% Allocate space for returned variables
S_output.x_t = zeros(d,S.nb_steps);

% Initial state
S_output.x_t(:,1) = S.initial(:,1);

% Generate the track
for t=2:S.nb_steps   
    % Propagate the state forward according to the given dynamic model
    S_output.x_t(:,t) = dynamic.model(S_output.x_t(:,t-1), dynamic);
end



