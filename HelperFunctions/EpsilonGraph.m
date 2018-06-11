function A = EpsilonGraph(X, epsilon)
%   Function to construct the epsilon-ball graph
%   For a given node, all other nodes within distance epsilon are its
%   neighbors
%
%   Input:
%       X: 2xN matrix of node positions
%       epsilon: scalar, maximum distance for neighbourhood
%
%   Output:
%       A: NxN adjacency matrix

N = size(X,1);
A = zeros(N,N);

% Compute distance between nodes
dist = distFast(X);

epsilon = epsilon*max(dist(:));

A(dist(:)<epsilon) = 1;
A(1:N+1:end) = 0;

