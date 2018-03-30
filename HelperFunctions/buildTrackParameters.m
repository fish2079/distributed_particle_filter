function [S, dynamic] = buildTrackParameters(sim_parameters)
%   Helper function to generate the necessary parameters that will be used
%   to construct the actual target track
%   We introduce this helper function so we can easily introduce new test
%   tracks 
%
%   Input:
%       sim_parameters: a struct that contains the basic parameters
%       such as track number
%
%   Output:
%       S: a struct that contains all the relevant
%       parameters for track construction such as number of sensors, area
%       length
%       dynamic: this struct contains all required parameters of target 
%       transition model

% First, we set the common parameters for all tracks
dynamic.model = @propagate_cv_ct; % dynamic model for target propagation
dynamic.T = 1; % sampling interval 
S.nb_steps = 50; % total number of time steps

% We use a big switch block set parameters according to the track number
% This allows us to easily introduce new track and modify existing ones
switch (sim_parameters.track)
    case 1
        S.nb_sensors = 4; % number of sensors 
        S.grid_sensor = true; % boolean flag to put the sensors in a grid
        S.initial = [45,45,5,0]';
        S.area_length = 75; % tracking area size = area_length^2 (km^2)
        S.area_xMin = 100;
        S.area_xMax = 175;
        S.area_yMin = 0;
        S.area_yMax = 75;
        S.broadcast_range = max(S.area_xMax-S.area_xMin, S.area_yMax-S.area_yMin)/(sqrt(S.nb_sensors)-1)*1.1;
        dynamic.a = 0.5; % turning rate
        dynamic.p = 0.05; % probability of turning
        dynamic.sigma_a = 0.05; % process noise standard deviation
    case 2
        S.nb_sensors = 16;
        S.grid_sensor = true;
        S.initial = [45,45,5,0]';
        S.area_length = 75;
        S.area_xMin = 50;
        S.area_xMax = 200;
        S.area_yMin = -50;
        S.area_yMax = 50;
        S.broadcast_range = max(S.area_xMax-S.area_xMin, S.area_yMax-S.area_yMin)/(sqrt(S.nb_sensors)-1)*1.1;
        dynamic.a = -0.25; 
        dynamic.p = 0.75;
        dynamic.sigma_a = 0.05;
    case 3
        S.nb_sensors = 9; 
        S.grid_sensor = true;
        S.initial = [45,45,5,0]';
        S.area_length = 75;
        S.area_xMin = S.initial(1)-S.area_length*0.5;
        S.area_xMax = S.initial(1)+S.area_length*0.5;
        S.area_yMin = S.initial(2)-S.area_length*0.5;
        S.area_yMax = S.initial(2)+S.area_length*0.5;
        S.broadcast_range = max(S.area_xMax-S.area_xMin, S.area_yMax-S.area_yMin)/(sqrt(S.nb_sensors)-1)*1.1;
        dynamic.a = 0.5; 
        dynamic.p = 0.95;
        dynamic.sigma_a = 0.05;
    case 4
        S.nb_sensors = 4; 
        S.grid_sensor = false;
        S.initial = [2,5,5,0]';
        S.area_length = 80;
        S.area_xMin = -40;
        S.area_xMax = 40;
        S.area_yMin = 20;
        S.area_yMax = 100;
        S.broadcast_range = 65; 
        dynamic.a = 0.5; 
        dynamic.p = 0.95;
        dynamic.sigma_a = 0.05;
end


% Track 3
% The initial target position is at the center of tracking area
% xMin = S.initial(1)-S.area_length*0.5;
% xMax = S.initial(1)+S.area_length*0.5;

% yMin = S.initial(2)-S.area_length*0.5;
% yMax = S.initial(2)+S.area_length*0.5;