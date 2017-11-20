function jointLikelihoodSums = likelihoodConsensus_real(z_n)

global parset;

switch parset.consensusAlgorithm(1)
    case 0
        jointLikelihoodSums = averageConsensus(z_n);
    case 1
        jointLikelihoodSums = quantizedCensi(z_n);
    case 2
        %TODO: update this function!
        %jointLikelihoodSums = quantizedAysal(z_n);
        disp('Aysal not implemented!!!');
    case 3
        %TODO: update this function!        
        %jointLikelihoodSums = pushSum(z_n,q_prime);
        disp('Push sum not implemented!!!');        
    otherwise
end

end