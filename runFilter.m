function [X_estimate, details] = runFilter(S, F, D, dynamic, obs)
%   Function to run one single trial using a single particle filter
%
%   Inputs:
%       S: Struct containing track parameters
%       F: Struct containing filter parameters
%       D: Struct containing measurement data
%       dynamics: Struct containing dynamic model parameters
%       obs: Struct containing measurement model paraleters
%
%   Output:
%       X_estimate: d x nb_steps matrix of estimated target states where d is
%       the target state dimension
%       position_error: 1 x nb_steps row vector of position estimation error
%       details: struct containing miscellaneous data for algorithm
%       debugging
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

% Allocate variables to store the final results
X_estimate = zeros(F.d, S.nb_steps); 

% Initialize particles
X_filterRepresentation(:,:,1) = F.initializeState(F);

% Initialize debug struct
details = struct();

% Iterate over each time step
for k = 1:S.nb_steps
    % Load data for this time step
    D_single.measurements = D.measurements{k};
    D_single.sensorID = D.sensorID{k};
    D_single.sensorLoc = D.sensorLoc{k};
        
    % Skip the predict step for time step 1
    if (k==1)
        F.initial = true;
    else
        F.initial = false;
    end
    
    % Run the particle filter for a single time step
    [X_filterRepresentation(:,:,k), details] = F.algorithm(X_filterRepresentation(:,:,max(k-1,1)), F, D_single, dynamic, obs, details);
    
    % State estimate (conditional expectation of particle distribution)
    X_estimate(:,k) = F.mmseStateEstimate(X_filterRepresentation(:,:,k));
    
    % Visualize the tracking results at individual time steps
    if (S.visualizeParticles && ~S.parallel)
        if (k==1)
            fig_handle = visualizeParticles(S, D_single, X_estimate(:,1:k), X_filterRepresentation(:,:,k), true, []);
        else
            fig_handle = visualizeParticles(S, D_single, X_estimate(:,1:k), X_filterRepresentation(:,:,k), true, fig_handle);
        end
        pause;
    end
end