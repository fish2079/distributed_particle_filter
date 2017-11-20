function [mu,sigma] = combinePartResults_avgConsensus(weightSum_part,mu_part,sigma_part)

global parset;

consensusNumOfIter = parset.consensusNumOfIter;

G   = parset.adjacencyMat;          %adjacenecy
D   = sum(G,2);                     %degrees
N   = size(G,1);

switch parset.consensusWeightModel(1)
    case 0 %constant step (uniform weights)
        eta = parset.consensusWeightModel(2);     %step size

        L = diag(D)-G;                            %Laplacian matrix
        P = eye(parset.numSensors)-eta*L;         %Perron matrix
    case 1 %maximum degree weights
        eta = 1/max(D);
        
        L = diag(D)-G;                            %Laplacian matrix
        P = eye(parset.numSensors)-eta*L;         %Perron matrix
        
    case 2 %local degree weights
        P = zeros(N);
        for i=1:N
            for j=(i+1):N
                if (G(i,j)==1)
                    P(i,j) = 1/(max(D(i),D(j)));
                end
            end
        end
        P = P+P';
        dg = diag(1-sum(P,2));
        P = P+dg;
    case 3 %metropolis weights
        P = zeros(N);
        for i=1:N
            for j=(i+1):N
                if (G(i,j)==1)
                    P(i,j) = 1/(1+max(D(i),D(j)));
                end
            end
        end
        P = P+P';
        dg = diag(1-sum(P,2));
        P = P+dg;
end


%Calculate the consensus values:

P_pow = P^consensusNumOfIter; %power of P matrix

weightSum = parset.numSensors*P_pow*weightSum_part.'; %sum the partial weight sums
mu = repmat(1./weightSum,1,4*parset.numTargets).*(parset.numSensors*P_pow*mu_part); %sum the partial means
%sum the partial covariances:
sigma =  zeros(size(sigma_part));
%do the summation for each component separately:
for r = 1:4*parset.numTargets
    for c = 1:4*parset.numTargets
        sigma(r,c,:) = reshape(1./weightSum.*(parset.numSensors*P_pow*reshape(sigma_part(r,c,:),parset.numSensors,1)),[1,1,parset.numSensors]);
    end
end

end