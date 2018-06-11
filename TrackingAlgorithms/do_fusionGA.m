function [MUave, Rave, scalars] = do_fusionGA(MUini,Rini_inv,F,...
    sensors_without_meas)

K = size(MUini,2);
% % nodes with measurements only
j = 1;
for i = 1:K
% Sensors without measurements will have MU and R of all zeros, we cannot
% have them in gossip biasing the averages. IF clause and also the variable
 % sensors_without_meas are taking care of this by assuring that only nodes
 % with measurements participate in gossip
    if nnz(Rini_inv(:,:,i))>0 && nnz(MUini(:,i))>0
        Rave_matrix(:,:,j) = Rini_inv(:,:,i);
        MUave_matrix(:,j) = Rini_inv(:,:,i) * MUini(:,i);
        j = j + 1;
    end
end

% [MUave_gossip,scalars1(1)] = gossip(MUave_matrix,F,sensors_without_meas);
% [MUave_max,scalars1(2)] = max_consensus_regular(MUave_gossip,F,sensors_without_meas);
MUave_max = sum(MUave_matrix,2);
MUave = MUave_max(:,1);

vectorized_Rave_matrix = reshape(Rave_matrix,F.d^2,K-length(sensors_without_meas));
% [Rave_gossip,scalars1(3)] = gossip(vectorized_Rave_matrix, F,sensors_without_meas);
% [Rave_max,scalars1(4)] = max_consensus_regular(Rave_gossip,F,sensors_without_meas);
Rave_max = sum(vectorized_Rave_matrix,2);
Rave = reshape(Rave_max,F.d,F.d);

% MUave = mean(MUave_matrix,2);
% Rave = mean(Rave_matrix,3);

Rave = chol_inv(Rave);

% min_eig = min(eig(Rave));
% if min_eig <= 1e-15
%     Rave =  Rave +  (-min_eig + 1e-15) * eye(size( Rave,2));
% end
MUave =  Rave *  MUave;

scalars = sum(scalars1);
%%% actually needs more transmission for determining #nodes with measurements !!!