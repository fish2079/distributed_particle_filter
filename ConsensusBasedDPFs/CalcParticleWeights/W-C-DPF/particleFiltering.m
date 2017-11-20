function x_n_hat = particleFiltering(z_n)
%Performs one PF step at each sensor. It uses 
%persistent variables to preserve the particles from the previous time 
%instant. As output, the function provides estimated targets positions and 
%velocities for each sensor.
%
%z_n - measurements for each sensor
%
%x_n_hat - numSensors x (4 * numTargets) matrix containing estimated
%position and velocity for each target and each sensor

global parset;

sqrtCu = sqrtm(parset.Cu); %square root of driving noise covariance (used for random number generation)

persistent particles_pos; %numParticles x (4*numTargets) matrix (contains positions of particles; the positions are the same for all sensors due to synchronized rnd generators)
%NOTE: We use persistent variables to keep the particles between calls for
%different time instants.

%Initialization:---------------------------------------------------------->

if isempty(particles_pos) %if empty, we are at the beginning of a new averaging cycle
    particles_pos = zeros(parset.numParticles,4*parset.numTargets);
    
    %sample from f(x0):
    for tt = 0:parset.numTargets-1
        particles_pos(:,tt*4+1:(tt+1)*4) = repmat((parset.x0(:,tt+1)).',parset.numParticles,1);
        particles_pos(:,tt*4+1:(tt+1)*4) = particles_pos(:,tt*4+1:(tt+1)*4)*parset.Phi.' + (randn(parset.numParticles,2)*parset.sqrtC_w0)*parset.Gamma0.';                    
    end
end

%<-------------------------------------------------------------------------

%Prediction step:--------------------------------------------------------->

%Sample from state transition pdf:
for tt = 0:parset.numTargets-1
    particles_pos(:,tt*4+1:(tt+1)*4) = particles_pos(:,tt*4+1:(tt+1)*4)*parset.Phi.' + (randn(parset.numParticles,2)*sqrtCu)*parset.Gamma.';
end

%<-------------------------------------------------------------------------

%Measurement update step:------------------------------------------------->
    
particles_weights_local = zeros(parset.numParticles,parset.numSensors);
for k = 1:parset.numSensors
    %singleLikelihood returns log of local likelihood...
    particles_weights_local(:,k) = singleLikelihood(particles_pos(:,repmat([1 1 0 0],1,parset.numTargets)==1),k,z_n(k));
end

switch parset.consensusType
    case 0 %perfect consensus
        particles_weights = sum(particles_weights_local,2);
    case 1 %realistic consensus
        particles_weights = parset.numSensors .* combineWeights_avgConsensus(particles_weights_local);
        
        %max/min consensuses so that all sensors obtain the same values:

        %we use the average of min and max weights...
        particles_weights = (max(particles_weights,[],2)+min(particles_weights,[],2))/2; %all sensors will use the same values -> this is need for resampling...
        
%         %we use the max weights...
%         particles_weights = max(particles_weights,[],2); %all sensors will use the same values -> this is need for resampling...
        
end

particles_weights = particles_weights - max(particles_weights);

particles_weights = exp(particles_weights);

%Normalization:
particles_weights = particles_weights ./ sum(particles_weights);

%Calculate target position estimate:
x_n_hat = repmat(sum(repmat(particles_weights,1,4*parset.numTargets).*particles_pos),parset.numSensors,1);

if parset.rougheningOnOff == 1 %roughening ON
    %Needed for roughening:
    E = max(particles_pos(:,1:4*parset.numTargets))-min(particles_pos(:,1:4*parset.numTargets)); %interval length between the maximum and minimum samples of each component (dimension) of the state vector
end        

%Resample:
particles_pos = particles_pos(systematicResamp(particles_weights,parset.numParticles,rand(1)),:);

if parset.rougheningOnOff == 1 %roughening ON    
    %Roughening:
    particles_pos(:,1:4*parset.numTargets) = roughening(particles_pos(:,1:4*parset.numTargets),E);
end
%<-------------------------------------------------------------------------

end