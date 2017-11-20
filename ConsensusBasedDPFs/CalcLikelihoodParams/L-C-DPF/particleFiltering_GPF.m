function x_n_hat = particleFiltering_GPF(predictionStep,jointLikelihoodSums)
%Performs one PF step (of the GPF algorithm) at each sensor with the 
%joint likelihood given by the sums obtained by likelihood consensus and 
%stored in jointLikelihoodSums. It uses persistent variables to preserve 
%the mean and covariance of the Gaussian posteriors from the previous time 
%instant. As output, the function provides estimated targets positions and
%velocities for each sensor.
%This version implements the Gaussian PF (from paper Kotecha, Djuric -
%Gaussian Particle Filtering, IEEE TSP, Oct. 2003)
%
%jointLikelihoodSums - matrix storing in rows the consensus sums defining 
%the joint likelihood for each sensor (or sensor measurements in the case 
%of no consensus (then it is a vector containing meas for all sensors); is 
%not processed only passed to jointLikelihood function
%predictionStep - if set to non-zero value it indicates that we only do
%particle prediction step; predicted particles are stored also in the
%global variable particlesApprox and are used to calculate reference points
%in Taylor approximation and for least squares approximation
%
%x_n_hat - numSensors x (4 * numTargets) matrix containing estimated
%position and velocity for each target and each sensor

global parset;
global particlesApprox; %used for approximation of measurement function

sqrtCu = chol(parset.Cu); %square root of driving noise covariance (used for random number generation)

%NOTE: We use persistent variable to keep the mean and covariance of the
%posterior between calls for different time instants.
persistent mu; %eq. (11)
persistent sigma; %eq. (11)

persistent particles_positions; %used to keep particles between calls of this function with predictionStep == 1 and predictionStep == 0

x_n_hat = zeros(parset.numSensors,4*parset.numTargets);

persistent targetPosIndx; %due to optimization

if (isempty(targetPosIndx))
    targetPosIndx = sparse(repmat([1 1 0 0],1,parset.numTargets));
end

if parset.consensusType ~= 0 %if LC is used (i.e. we do not simulate centralized PF)
%Only in the case of LC we do PF for each sensor. In the case of no LC 
%(centraziled PF) we do processing only for 1 sensor (fusion center) with
%an exact joint likelihood that contains the measurements of all sensors.

if predictionStep == 1 %do prediction step

%Initialization:---------------------------------------------------------->
if isempty(mu) %if empty, we are at the beginning of a new averaging cycle
    %We sample particles from the prior since there is no previous Gaussian
    %posterior available.
    particles_positions = zeros(parset.numParticles,4*parset.numTargets,parset.numSensors);
    
    %sample from f(x0):
    for tt = 0:parset.numTargets-1
        particles_positions(:,tt*4+1:(tt+1)*4,:) = repmat((parset.x0(:,tt+1)).',[parset.numParticles,1,parset.numSensors]);
        for k = 1:parset.numSensors
%             particles_positions(:,tt*4+1:(tt+1)*4,k) = particles_positions(:,tt*4+1:(tt+1)*4,k) + randn(parset.numParticles,4)*parset.sqrtC0;
            particles_positions(:,tt*4+1:(tt+1)*4,k) = particles_positions(:,tt*4+1:(tt+1)*4,k)*parset.Phi.' + (randn(parset.numParticles,2)*parset.sqrtC_w0)*parset.Gamma0.';            
        end
    end    
%<-------------------------------------------------------------------------
else %sample particles from the previous Gaussian posterior
        
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

for k = 1:parset.numSensors
    %Sample from state transition pdf:
    for tt = 0:parset.numTargets-1
        particles_positions(:,tt*4+1:(tt+1)*4,k) = particles_positions(:,tt*4+1:(tt+1)*4,k)*parset.Phi.' + (randn(parset.numParticles,2)*sqrtCu)*(parset.Gamma.');    
    end
end

%NOTE: The predicted particles are stored in the presistent variable. We 
%therefore don't need to repeat the prediction step when this function is 
%called again with predictionStep = 0.

particles =  1/parset.numParticles*ones(parset.numParticles,4*parset.numTargets+1,parset.numSensors); %we need to store in the particlesApprox also the weights...
particles(:,1:4*parset.numTargets,:) = particles_positions;
particlesApprox = particles; %used for approximation of measurement function

else %predictionStep == 0 -- do measurement update step
%Do PF for each sensor:--------------------------------------------------->

%NOTE: particles_positions already contains the predicted particles (we did
%the prediction step in the previous call with predictionStep == 1).

mu = zeros(parset.numSensors,4*parset.numTargets);
sigma = zeros(4*parset.numTargets,4*parset.numTargets,parset.numSensors);

if parset.reducedComplex == 1 %with complexity reduction
    weightSum_partial = zeros(1,parset.numSensors); %used to store the partial sum of weights for each sensor
end

for k = 1:parset.numSensors

    % NOTE: This function is always called using predictionStep = 1 prior to
    % calling it with predictionStep = 0 (see simulate.m). The persistent 
    % variable "particles_positions" already contains the predicted 
    % particles  => there is no need to do the prediction step here.
    
    %Calculate particle weights:
    particles_k_weights = jointLikelihood(particles_positions(:,targetPosIndx==1,k),jointLikelihoodSums(k,:));
    
    if parset.reducedComplex == 1 %with complexity reduction
        %If we use complexity reduction, the weights are not normalized, we
        %calculate their partial sum and we compute the not-normalized
        %means and covariances.
        weightSum_partial(k) = sum(particles_k_weights);
    end
    
    %Calculate variance of particles:
    tmp1 = repmat(particles_k_weights,1,4*parset.numTargets).*particles_positions(:,:,k);
    sigma(:,:,k) =  tmp1'*particles_positions(:,:,k);
    
    %Calculate the mean of particles:
    mu(k,:) = sum(tmp1);
    
    if parset.reducedComplex == 0 %no complexity reduction
        %do weight normalization:
        weightSum = sum(particles_k_weights);
        mu(k,:) = 1/weightSum * mu(k,:);
        sigma(:,:,k) = 1/weightSum * sigma(:,:,k);
    end
    
end

%Combination of results:------------------------------------------>
%In case we simulate with "the reduced complexity", we need to
%combine the results (means & covariances) from individual sensors in order
%to obtain the full-blown estimate.
if parset.reducedComplex == 1 %with complexity reduction
    if parset.consensusType == 1 %perfect consensus
        weightSum_fullBlown = sum(weightSum_partial);
        mu = repmat(1/weightSum_fullBlown * sum(mu),parset.numSensors,1);
        sigma = repmat(1/weightSum_fullBlown * sum(sigma,3),[1,1,parset.numSensors]);       
    elseif parset.consensusType == 2 %realistic consensus
        %average consensus is used:
        [mu,sigma] = combinePartResults_avgConsensus(weightSum_partial,mu,sigma); 
    end
end
%<-----------------------------------------------------------------

%subtract the "mean" to obtain covariance:
for k = 1:parset.numSensors 
    sigma(:,:,k) = sigma(:,:,k) - mu(k,:).'*mu(k,:);
end

%Target position estimate:
x_n_hat = mu;

%<-------------------------------------------------------------------------

end

%==========================================================================

else %parset.consensusType == 0; no LC (centralized PF) -> we only need to do 1 PF (at fusion center)
    
if predictionStep == 1 %do prediction step
    
%Initialization:---------------------------------------------------------->
if isempty(mu) %if empty, we are at the beginning of a new averaging cycle
%We sample particles from the prior since there is no previous Gaussian
%posterior available.
    particles_positions = zeros(parset.numParticles,4*parset.numTargets);
      
    %sample from f(x0):
    for tt = 0:parset.numTargets-1
        particles_positions(:,tt*4+1:(tt+1)*4) = repmat((parset.x0(:,tt+1)).',parset.numParticles,1);
%         particles_positions(:,tt*4+1:(tt+1)*4) = particles_positions(:,tt*4+1:(tt+1)*4) + randn(parset.numParticles,4)*parset.sqrtC0;
        particles_positions(:,tt*4+1:(tt+1)*4) = particles_positions(:,tt*4+1:(tt+1)*4)*parset.Phi.' + (randn(parset.numParticles,2)*parset.sqrtC_w0)*parset.Gamma0.';
    end    
%<-------------------------------------------------------------------------
else %sample particles from the previous Gaussian posterior
        
%     sqrtSig = chol(sigma);
%     particles_positions = repmat(mu,parset.numParticles,1) + randn(parset.numParticles,4*parset.numTargets)*sqrtSig;                
    try
        sqrtSig = chol(sigma);
        particles_positions = repmat(mu,parset.numParticles,1) + randn(parset.numParticles,4*parset.numTargets)*sqrtSig;                    
    catch %#ok<CTCH>
%             disp('WARNING in particleFiltering: Diagonalizing singular covariance matrix!');
        try
            sqrtSig = chol(diag(diag(sigma)));
            particles_positions = repmat(mu,parset.numParticles,1) + randn(parset.numParticles,4*parset.numTargets)*sqrtSig;
        catch %#ok<CTCH>
%                 disp('WARNING in particleFiltering: Singular covariance matrix!');
            %If the covariance is singular and close to zero, it means that the Gaussian is
            %almost a Dirac -> all samples drawn will be the same -> we set
            %them to be equal to the mean.
            particles_positions = repmat(mu,parset.numParticles,1);
        end
    end

end    

%Prediction step:

%Sample from state transition pdf:
for tt = 0:parset.numTargets-1
    particles_positions(:,tt*4+1:(tt+1)*4) = particles_positions(:,tt*4+1:(tt+1)*4)*parset.Phi.' + (randn(parset.numParticles,2)*sqrtCu)*parset.Gamma.';
end

else %predictionStep == 0 -- do measurement update step

%NOTE: particles_positions already contains the predicted particles (we did
%the prediction step in the previous call with predictionStep == 1).
            
%Calculate particle weights:
particles_weights = jointLikelihood(particles_positions(:,targetPosIndx==1),jointLikelihoodSums);
%NOTE: jointLikelihoodSums now contains measurements of all sensors

%Calculate the mean of particles:
mu = 1/sum(particles_weights) * sum(repmat(particles_weights,1,4*parset.numTargets).*particles_positions);

%Calculate variance of particles:
sigma = zeros(4*parset.numTargets);
for p = 1:parset.numParticles
    sigma = sigma + particles_weights(p)*particles_positions(p,:).'* particles_positions(p,:);
end
sigma = 1/sum(particles_weights) * sigma - mu.'*mu;

%Target position estimate:
x_n_hat = repmat(mu,parset.numSensors,1);

end    
end

end