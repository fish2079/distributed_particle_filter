function newParticles = roughening(oldParticles)
%This function performs roughening on the set of old (resampled) particles 
%to prevent sample empoverishment. It adds a random jitter (Gaussian noise)
%to each particle.
%See Gordon et al. - Novel approach to nonlinear/non-Gaussian Bayesian
%state estimation (pg. 112 botom left)
%
%oldParticles - set of old particles; matrix containing in rows the
%particles (4*numTargets cols -> no weights)
%
%newParticles - new set of particles (contains "jittered" old particles)

global parset;

numParticles = size(oldParticles,1);
newParticles = zeros(numParticles,4*parset.numTargets);

E = max(oldParticles(:,1:4*parset.numTargets))-min(oldParticles(:,1:4*parset.numTargets)); %interval length between the maximum and minimum samples of each component (dimension) of the state vector


for tt = 0:parset.numTargets-1
%check if interval lengths are non-zero:
if E(tt*4+3) == 0 && E(tt*4+4) == 0
    %NOTE: If all interval lengths are zero, it means that we only have 
    %one particle (or multiple copies of one particle). Calculation of 
    %sigma in the standard way leads to zero sigma values. This in turn
    %causes that no jitter is added to the multiple copies of one particle
    %and hence the output of this function is the same as its input.
    %Therefore we set sigma to the value of K. This way, roughening returns
    %different particles also in the case when the input contains only 
    %multiple copies of a single particle. This in turn helps in clustering
    %functions.
    E(tt*4+3) = 10^-2;
    E(tt*4+4) = 10^-2;
end
end

sigma = parset.rougheningTuningK*E*numParticles^(-1/4);%standard deviation of the jitter noise for each component direction

for tt = 0:parset.numTargets-1
     newParticles(:,tt*4+1:(tt+1)*4) = oldParticles(:,tt*4+1:(tt+1)*4) + (randn(parset.numParticles,2)*diag([sigma(4*tt+3),sigma(4*tt+4)]))*parset.Gamma.';
end

end