function position_error = computePositionError(x_est, x_true)
%   Function to compute the position estimation error. The error is
%   computed as L2-norm of discrepancy between estimate and true position
%
%   Input:
%       x_est: dxK matrix of estimated states, d is the target state
%       dimension and the 2 first entries denote the positions, K is the
%       total numer of estimates/time steps/particles
%       x_true: dxK matrix of true states
%
%   Output: 
%       position_error: 1xK row vector of position estimation error
%

dif = x_est-x_true;

dif_position = dif(1:2,:);

dif_pos_squared = dif_position.^2;

position_error = sqrt(sum(dif_pos_squared,1));