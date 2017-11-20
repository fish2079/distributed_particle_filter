function x_n_hat = particleFiltering_SIR(predictionStep,jointLikelihoodSums)
%Performs one PF step (of the SIR PF algorithm) at each sensor with the 
%joint likelihood given by the sums obtained by likelihood consensus and 
%stored in jointLikelihoodSums. It uses persistent variables to preserve 
%the particles from the previous time instant. As output, the function 
%provides estimated targets positions and velocities for each sensor.
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

sqrtCu = sqrtm(parset.Cu); %square root of driving noise covariance (used for random number generation)

persistent particles; %numParticles x (4*numTargets + 1) x numSensors matrix (first 4*numTargets columns store particles; the last one stores weights)
%NOTE: We use persistent variable to keep the particles between calls for
%different time instants.
persistent targetPosIndx; %due to optimization

x_n_hat = zeros(parset.numSensors,4*parset.numTargets);

if (isempty(targetPosIndx))
    targetPosIndx = sparse(repmat([1 1 0 0],1,parset.numTargets));
end

if parset.consensusType ~= 0 %if LC is used (i.e. we do not simulate centralized PF)
%Only in the case of LC we do PF for each sensor. In the case of no LC 
%(centraziled PF) we do processing only for 1 sensor (fusion center) with
%an exact joint likelihood that contains the measurements of all sensors.

%Initialization:---------------------------------------------------------->
if isempty(particles) %if empty, we are at the beginning of a new averaging cycle
    particles = zeros(parset.numParticles,4*parset.numTargets+1,parset.numSensors);
    
    %sample from f(x0):
    for tt = 0:parset.numTargets-1
        particles(:,tt*4+1:(tt+1)*4,:) = repmat((parset.x0(:,tt+1)).',[parset.numParticles,1,parset.numSensors]);
        for k = 1:parset.numSensors
%             particles(:,tt*4+1:(tt+1)*4,k) = particles(:,tt*4+1:(tt+1)*4,k) + randn(parset.numParticles,4)*parset.sqrtC0;
            particles(:,tt*4+1:(tt+1)*4,k) = particles(:,tt*4+1:(tt+1)*4,k)*parset.Phi.' + (randn(parset.numParticles,2)*parset.sqrtC_w0)*parset.Gamma0.';            
        end
    end
    
    %Set the weights of the particles:
    particles(:,4*parset.numTargets+1,:) = 1/parset.numParticles;
end
%<-------------------------------------------------------------------------

if predictionStep == 1 %do prediction step

for k = 1:parset.numSensors
    %Sample from state transition pdf:
    for tt = 0:parset.numTargets-1
        particles(:,tt*4+1:(tt+1)*4,k) = particles(:,tt*4+1:(tt+1)*4,k)*parset.Phi.' + (randn(parset.numParticles,2)*sqrtCu)*(parset.Gamma.');    
    end
end

%NOTE: The predicted particles are stored in the presistent variable. We 
%therefore don't need to repeat the prediction step when this function is 
%called again with predictionStep = 0.

particlesApprox = particles; %used for approximation of measurement function

else %predictionStep == 0 -- do measurement update step
    
%Do PF for each sensor:--------------------------------------------------->

for k = 1:parset.numSensors

% NOTE: This function is always called using predictionStep = 1 prior to
% calling it with predictionStep = 0 (see simulate.m). The persistent 
% variable "particles" already contains the predicted particles (with old
% weights) => there is no need to do the prediction step here.

    %Calculate particle weights:
    particles(:,4*parset.numTargets+1,k) = particles(:,4*parset.numTargets+1,k) .* jointLikelihood(particles(:,targetPosIndx==1,k),jointLikelihoodSums(k,:));
    particles(:,4*parset.numTargets+1,k) = particles(:,4*parset.numTargets+1,k)/sum(particles(:,4*parset.numTargets+1,k));
   
    %Calculate target position estimate:
    x_n_hat(k,:) = sum(repmat(particles(:,4*parset.numTargets+1,k),1,4*parset.numTargets).*particles(:,1:4*parset.numTargets,k));
    
    %Resample
    particles(:,1:4*parset.numTargets,k) = particles(systematicResamp(particles(:,4*parset.numTargets+1,k),parset.numParticles,rand(1)),1:4*parset.numTargets,k);
    particles(:,4*parset.numTargets+1,k) = 1/parset.numParticles;
    
    if parset.rougheningOnOff == 1 %roughening ON    
        %Roughening:
        particles(:,1:4*parset.numTargets,k) = roughening(particles(:,1:4*parset.numTargets,k));
    end
    
end

%<-------------------------------------------------------------------------

end

%==========================================================================

else %parset.consensusType == 0; no LC (centralized PF) -> we only need to do 1 PF (at fusion center)
    
%Initialization:---------------------------------------------------------->
if isempty(particles) %if empty, we are at the beginning of a new averaging cycle
    particles = zeros(parset.numParticles,4*parset.numTargets+1);
    
    %sample from f(x0):
    for tt = 0:parset.numTargets-1
        particles(:,tt*4+1:(tt+1)*4) = repmat((parset.x0(:,tt+1)).',parset.numParticles,1);
%         particles(:,tt*4+1:(tt+1)*4) = particles(:,tt*4+1:(tt+1)*4) + randn(parset.numParticles,4)*parset.sqrtC0;
        particles(:,tt*4+1:(tt+1)*4) = particles(:,tt*4+1:(tt+1)*4)*parset.Phi.' + (randn(parset.numParticles,2)*parset.sqrtC_w0)*parset.Gamma0.';
    end
    
    %Set the weights of the particles:
    particles(:,4*parset.numTargets+1) = 1/parset.numParticles;
end
%<-------------------------------------------------------------------------

if predictionStep == 1 %do prediction step

    %Sample from state transition pdf:
    for tt = 0:parset.numTargets-1
        particles(:,tt*4+1:(tt+1)*4) = particles(:,tt*4+1:(tt+1)*4)*parset.Phi.' + (randn(parset.numParticles,2)*sqrtCu)*parset.Gamma.';
    end

else %predictionStep == 0 -- do measurement update step

%NOTE: the variable particles already contains the predicted particles (we 
%did the prediction step in the previous call with predictionStep == 1).
    
    %Calculate particle weights:
    particles(:,4*parset.numTargets+1) = particles(:,4*parset.numTargets+1) .* jointLikelihood(particles(:,targetPosIndx==1),jointLikelihoodSums);
    particles(:,4*parset.numTargets+1) = particles(:,4*parset.numTargets+1)/sum(particles(:,4*parset.numTargets+1));
   
    %Calculate target position estimate:
    x_n_hat = repmat(sum(repmat(particles(:,4*parset.numTargets+1),1,4*parset.numTargets).*particles(:,1:4*parset.numTargets)),parset.numSensors,1);
    
    %Resample
    particles(:,1:4*parset.numTargets) = particles(systematicResamp(particles(:,4*parset.numTargets+1),parset.numParticles,rand(1)),1:4*parset.numTargets);
    particles(:,4*parset.numTargets+1) = 1/parset.numParticles;
    
    if parset.rougheningOnOff == 1 %roughening ON    
        %Roughening:
        particles(:,1:4*parset.numTargets) = roughening(particles(:,1:4*parset.numTargets));
    end

end
    
end

end