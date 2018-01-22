function [x_aggregate, aggregate_error_ratio, initial_error_ratio] = computeAggregateGossip(x_initial, A, max_gossip_iter)
%   Function to compute the aggregate values using gossip
%   A consensus round is run in the end to ensure all sensors have the same
%   values
%
%   Input:
%       x_initial: column vector of initial values. If x_initial is a
%       matrix, each column represents one set of values to be aggregated
%       A: NxN adjacency matrix where N is the number of sensors
%       max_gossip_iter: scalar, maximum number of gossip iterations
%   
%   Output:
%       x_aggregate: column vector of sensor values after distributed
%       summation
%       aggregate_error: row vector of aggregate error due to using gossip
%       and max consensus

% First, compute the Metropolis weight for the aggregate
% Compute the degree of each sensor
D = sum(A,1);
N=size(A,1);
W = zeros(N,N);
for i=1:N
    for j=1:N
        if (i~=j && A(i,j)==1)
            W(i,j) = 1/(1+max(D(i),D(j)));
        end
    end
end

W(1:N+1:end) = 1-sum(W,2);

W_final = W^max_gossip_iter;

x_aggregate_gossip = W_final*(N*x_initial);

x_aggregate = max(x_aggregate_gossip,[],1);

% x_aggregate_true = sum(x_initial,1);
% aggregate_error_ratio = (x_aggregate - x_aggregate_true)./sqrt(sum((N*x_initial-x_aggregate_true).^2,1));
% initial_error_ratio = abs(sqrt(sum((N*x_initial-x_aggregate_true).^2,1))./sum(x_initial,1));

% aggregate_error_ratio = [aggregate_error_ratio; initial_error_ratio];
aggregate_error_ratio = abs((x_aggregate(1,:)-sum(x_initial,1))./sum(x_initial,1));
yo=5;


