warning('off','all');
clear;clc;

path = 'Results_N\Track1\';
% Number of particles for the filter
N_vector = [100, 250, 500, 750, 1000];

% Number of random trials
sim_parameters.no_trials = 500; 

RMSE = [];
RMSEFull = [];
runtime = [];
runtimeFull = [];

% Loop through each choice of particle number
for i=1:numel(N_vector)
    % Set number of particles
    sim_parameters.N = N_vector(i); 
    
    % Store the tracking results
    filename{i} = [path, 'Track1_N',num2str(sim_parameters.N),'_trials',num2str(sim_parameters.no_trials),'.mat'];
    
    [RMSE_sf, RMSEFull_sf, runtime_sf, runtimeFull_sf, detail{i}] = extractResults(filename{i});
    
    RMSE = [RMSE; RMSE_sf];
    RMSEFull = [RMSEFull; RMSEFull_sf];
    runtime = [runtime; runtime_sf];
    runtimeFull = [runtimeFull; runtimeFull_sf];
end

