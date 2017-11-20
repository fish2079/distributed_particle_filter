function lh = singleLikelihood(x_n,sensor,z)
%Computes likelihood values for each input sample for the given sensor and
%its measurement.
%
%x_n - numSamples x 2*numTargets vector containing the samples where the likelihood
%is to be evaluated
%sensor - index of sensor for which likelihood is to be evaluated
%z - measurement of the sensor
%
%lh - numSamples x 1 vector of likelihood values evaluated at x_n
%(likelihood is evaluated up to a multiplicative constant factor)

global parset;

sensorPos = parset.sensorsPos(sensor,:);%to speed up execution

switch parset.measModel
    case 0 %inverse decay
        sx = size(x_n,1);%number of samples

        A_r_cumulative = zeros(sx,1);
        for tt = 0:parset.numTargets-1
            r = x_n(:,tt*2+1:(tt+1)*2) - repmat(sensorPos,sx,1);
            r = hypot(r(:,1),r(:,2));
            r(r<parset.d0) = parset.d0; %due to numerical stability... (if the distance is smaller than d0, we set it to d0; this is equivalent to setting A/r equal to A/d0)
            A_r_cumulative = A_r_cumulative + parset.A./(r.^parset.invPow);
        end        
        lh = -((z - A_r_cumulative).^2 )* 0.5 / (parset.w_sigma^2);
        
        lh_old = lh; %due to numerical instability for certain settings of the model 
        lh = lh.*parset.numSensors;
        
%     case 1 %exponential decay
%         %TODO
%     case 2 %bearing-only
%         lh = -(0.5/(parset.w_sigma^2))*(z - atan2(x_n(:,2)-sensorPos(2),x_n(:,1)-sensorPos(1))).^2;    
end

lh = lh - max(lh);
lh = exp(lh);

%due to numerical instability for certain settings of the model:---------->

N_eff = 1/sum((lh./sum(lh)).^2); %determine "effective" number of particles (number of particles with signifficant weights...)
    
if N_eff < 0.001 * parset.numParticles
    
    %disp('Too few effective particles...')
    lh = lh_old;
   
    lh = lh - max(lh);
    lh = exp(lh);
end
%<-------------------------------------------------------------------------

end