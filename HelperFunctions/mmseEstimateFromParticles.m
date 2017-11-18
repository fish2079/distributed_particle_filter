function x_estimate = mmseEstimateFromParticles(X_filterRepresentation)
%   Compute the minimum mean-squared-error state estimate from particle cloud
%
%   Inputs:
%       X_filterRepresentation:  (d+1)xN particle cloud matrix 
%
%   Outputs:
%       x_estimate:  d-by-1 state estimate vector
%
% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

d = size(X_filterRepresentation,1)-1;
x_estimate = X_filterRepresentation(1:d,:)*X_filterRepresentation(d+1,:)';

