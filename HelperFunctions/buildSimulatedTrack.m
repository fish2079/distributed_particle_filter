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

%% Place sensors in the tracking area
xMin = S.area_xMin;
xMax = S.area_xMax;
yMin = S.area_yMin;
yMax = S.area_yMax;
area_x_length = xMax - xMin;
area_y_length = yMax - yMin;
% Place sensors in a grid if needed
if (S.grid_sensor)
    rootK = sqrt(S.nb_sensors);
    sensorX = xMin:(xMax-xMin)/(rootK-1):xMax;
    sensorY = yMin:(yMax-yMin)/(rootK-1):yMax;   
    [sensorX, sensorY] = meshgrid(sensorX, sensorY);
    sensorPos = [sensorX(:), sensorY(:)]';
else
    % Alternatively, place sensors randomly in the tracking field
    offset = rand(S.nb_sensors, 2);
    sensorPos = [xMin + offset(:,1)*area_x_length, yMin + offset(:,2)*area_y_length]';
end
% Store all sensor positions
S_output.sensorPos = sensorPos;

% Construct sensor communication adjacency matrix
A = zeros(S.nb_sensors);
range = S.broadcast_range;
for i=1:S.nb_sensors
    for j=i+1:S.nb_sensors
        A(i,j) = (norm(sensorPos(:,i)-sensorPos(:,j))<range);
        A(j,i) = A(i,j);
    end
end
S_output.A = A;

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



