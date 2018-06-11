function [inverseM] = chol_inv(M)

symmetricM = mean(cat(3,M,M'),3);
[cholM,p] = chol(symmetricM);
if p>0
    warning('matrix is not positive definite')
end
cholInvM =  cholM \ eye(size(cholM));
inverseM = cholInvM*cholInvM';