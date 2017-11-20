function x_n_hat = particleFiltering(z_n)
%x_n_hat - numSensors x 4 matrix containing estimated target position for 
%each sensor in rows

global parset;

sqrtCu = chol(parset.Cu); %square root of driving noise covariance (used for random number generation)

%NOTE: We use persistent variable to keep the mean and covariance of the
%posterior between calls for different time instants.
persistent mu;
persistent sigma;

persistent targetPosIndx; %due to optimization

if (isempty(targetPosIndx))
    targetPosIndx = sparse(repmat([1 1 0 0],1,parset.numTargets));
end

%Initialization:---------------------------------------------------------->
if isempty(mu) %if empty, we are at the beginning of a new averaging cycle

    %We sample particles from the prior since there is no previous Gaussian
    %posterior available.
    particles_positions = zeros(parset.numParticles,4*parset.numTargets,parset.numSensors);
    
    %sample from f(x0):
    for tt = 0:parset.numTargets-1
        particles_positions(:,tt*4+1:(tt+1)*4,:) = repmat((parset.x0(:,tt+1)).',[parset.numParticles,1,parset.numSensors]);
        for k = 1:parset.numSensors
            particles_positions(:,tt*4+1:(tt+1)*4,k) = particles_positions(:,tt*4+1:(tt+1)*4,k)*parset.Phi.' + (randn(parset.numParticles,2)*parset.sqrtC_w0)*parset.Gamma0.';            
        end
    end  
%<-------------------------------------------------------------------------
else
    %sample particles from the Gaussian posterior
    particles_positions = zeros(parset.numParticles,4*parset.numTargets,parset.numSensors); %we don't need the weights column
    
    for k = 1:parset.numSensors
%         sqrtSig = chol(sigma(:,:,k));        
%         particles_positions(:,:,k) = repmat(mu(k,:),parset.numParticles,1) + randn(parset.numParticles,4*parset.numTargets)*sqrtSig;                
        try
            sqrtSig = chol(sigma(:,:,k));
            particles_positions(:,:,k) = repmat(mu(k,:),parset.numParticles,1) + randn(parset.numParticles,4*parset.numTargets)*sqrtSig;        
        catch %#ok<CTCH>
%             disp('WARNING in particleFiltering: Diagonalizing singular covariance matrix!');
            try
                sqrtSig = chol(diag(diag(sigma(:,:,k))));
                particles_positions(:,:,k) = repmat(mu(k,:),parset.numParticles,1) + randn(parset.numParticles,4*parset.numTargets)*sqrtSig;                
            catch %#ok<CTCH>
%                 disp('WARNING in particleFiltering: Singular covariance matrix!');
                %If the covariance is singular and close to zero, it means that the Gaussian is
                %almost a Dirac -> all samples drawn will be the same -> we set
                %them to be equal to the mean.
                particles_positions(:,:,k) = repmat(mu(k,:),parset.numParticles,1);
            end
        end
        %TODO: Throw an exception instead of "correcting" the covariance
        %matrix? This might be reasonable for settings where the case of
        %singular covariance is very rare.
    end     
end

%Do auxiliary PF for each sensor:----------------------------------------->

% %AUX PF:
% muk = particles_positions; %the mean of the state transition pdf
% lambda = zeros(parset.numSensors,parset.numParticles);
% for k = 1:parset.numSensors 
%     %Sample from state transition pdf:
%     for tt = 0:parset.numTargets-1
%         muk(:,tt*4+1:(tt+1)*4,k) = muk(:,tt*4+1:(tt+1)*4,k)*parset.Phi.';    
%     end
%     lambda(k,:) = singleLikelihood(muk(:,targetPosIndx==1,k),k,z_n(k));
%     lambda(k,:) = lambda(k,:)/sum(lambda(k,:));
%     indices = systematicResamp(lambda(k,:),parset.numParticles,rand(1));
%     particles_positions(:,1:4*parset.numTargets,k) = particles_positions(indices,1:4*parset.numTargets,k);
%     lambda(k,:) = lambda(k,indices);
% end

for k = 1:parset.numSensors 
    %Sample from state transition pdf:
    for tt = 0:parset.numTargets-1
        particles_positions(:,tt*4+1:(tt+1)*4,k) = particles_positions(:,tt*4+1:(tt+1)*4,k)*parset.Phi.' + (randn(parset.numParticles,2)*sqrtCu)*(parset.Gamma.');    
    end  
end

mu = zeros(parset.numSensors,4*parset.numTargets);
sigma = zeros(4*parset.numTargets,4*parset.numTargets,parset.numSensors);

for k = 1:parset.numSensors
           
    %Calculate particle weights:
    particles_k_weights = singleLikelihood(particles_positions(:,targetPosIndx==1,k),k,z_n(k));
%     particles_k_weights = particles_k_weights./(lambda(k,:)).';%AUX PF
    particles_k_weights = particles_k_weights./sum(particles_k_weights);
    
    %Calculate mean and variance of particles:
    tmp1 = repmat(particles_k_weights,1,4*parset.numTargets).*particles_positions(:,:,k);
    mu(k,:) = sum(tmp1);
    mu_rep = repmat(mu(k,:),parset.numParticles,1);
    tmp1 = repmat(particles_k_weights,1,4*parset.numTargets).*(particles_positions(:,:,k)-mu_rep);
    sigma(:,:,k) = tmp1'*(particles_positions(:,:,k)-mu_rep);
    
end

%Combination of local results:-------------------------------------------->
invSigma = zeros(size(sigma));
mu_tmp = zeros(size(mu));
for k = 1:parset.numSensors
    lastwarn('');    
    invSigma(:,:,k) = inv(sigma(:,:,k));
    if ~isempty(lastwarn) %if warnings are encountered during matrix inversion, we repeat the averaging cycle
        error(lastwarn);
        lastwarn('');
    end   
    mu_tmp(k,:) = (invSigma(:,:,k)*mu(k,:).').';
end

if parset.consensusType == 0 %perfect consensus
    lastwarn('');    
    sigma = repmat(inv((1/parset.numSensors).*sum(invSigma,3)),[1,1,parset.numSensors]);
    if ~isempty(lastwarn) %if warnings are encountered during matrix inversion, we repeat the averaging cycle
        error(lastwarn);
        lastwarn('');
    end    
    mu = repmat(((1/parset.numSensors).*sum(mu_tmp))*sigma(:,:,1),parset.numSensors,1);
elseif parset.consensusType == 1 %realistic consensus
    %average consensus is used:
    [mu_tmp,invSigma] = combinePartResults_avgConsensus(mu_tmp,invSigma);
    for k = 1:parset.numSensors
        lastwarn('');        
        sigma(:,:,k) = inv(invSigma(:,:,k));
        if ~isempty(lastwarn) %if warnings are encountered during matrix inversion, we repeat the averaging cycle
            error(lastwarn);
            lastwarn('');
        end        
        mu(k,:) = (sigma(:,:,k)*mu_tmp(k,:).').';
    end
end
%<-------------------------------------------------------------------------

%Target position estimate:
x_n_hat = mu;

%<-------------------------------------------------------------------------

end