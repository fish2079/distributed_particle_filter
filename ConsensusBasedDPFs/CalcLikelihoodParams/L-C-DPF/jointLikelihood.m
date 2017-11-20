function lh = jointLikelihood(x_n,jls)

%Computes joint likelihood values for each input sample.
%
%x_n - numSamples x 2*numTargets vector containing the samples where the 
%likelihood is to be evaluated
%jls - vector containing the "joint likelihood sums" computed via consensus
%that completely describes the joint likelihood for given sensor; for the 
%no-consensus case this vector contains measurement for each sensor (the 
%meas. are used to calculate true joint likelihood in a centralized way)
%
%lh - numSamples x 1 vector of likelihood values evaluated at x_n
%(likelihood is evaluated up to a multiplicative constant factor)

global parset;

if parset.consensusType == 0
%No consensus:------------------------------------------------------------>
%True joint likelihood for purposes of centralized PF is calculated

%NOTE: jls contains measurement for each sensor
sx = size(x_n,1); %number of samples

numSensors = parset.numSensors; %to speed up execution
sensorsPos = parset.sensorsPos; %to speed up execution

A = parset.A; %to speed up execution
invPow = parset.invPow; %to speed up execution
w_sigma = parset.w_sigma; %to speed up execution

switch parset.measModel
    case 0 %inverse decay
        A_r_cumulative = zeros(sx*parset.numSensors,1); %due to the sum over all targets in measurement function...
        if parset.trackLossDetOnOff > 1
            A_r_cumulative_detection = zeros(sx*parset.numSensors,parset.numTargets);
        end        
        for tt = 0:parset.numTargets-1
            r = reshape(repmat(reshape(x_n(:,tt*2+1:(tt+1)*2),1,2*sx),numSensors,1),sx*numSensors,2)-repmat(sensorsPos,sx,1);
            r = hypot(r(:,1),r(:,2));
            r(r<parset.d0) = parset.d0; %due to numerical stability... (if the distance is smaller than d0, we set it to d0; this is equivalent to setting A/r equal to A/d0)
            A_r_cumulative = A_r_cumulative + A./(r.^invPow);
            if parset.trackLossDetOnOff > 1
                A_r_cumulative_detection(:,tt+1) = A./(r.^invPow);
            end            
        end
        
        if parset.trackLossDetOnOff == 2
            %diregard sensors that are too close to the target...
            
            A_r_cumulative = reshape(A_r_cumulative,parset.numSensors,sx);
            
            droppedMeas = zeros(parset.numSensors,1) == 1; %vector of logicals initialized to 0
            for tt = 0:parset.numTargets-1
                A_r_cumulative_detection_tmp = reshape(A_r_cumulative_detection(:,tt+1),parset.numSensors,sx);
                droppedMeas = droppedMeas + any(A_r_cumulative_detection_tmp > (1/(parset.trackLossDetTrshld^parset.invPow))*parset.A,2);
            end      
            droppedMeas = droppedMeas > 0; %convert back to logicals...
            
            A_r_cumulative = A_r_cumulative(~droppedMeas.',:);
            
            numSensors = sum(~droppedMeas);
            
            A_r_cumulative = reshape(A_r_cumulative,numSensors*sx,1);
            
            jls = jls(~droppedMeas);
                        
%             if sum(droppedMeas) > 0
%                 disp('DROPPING');
%             end
        end
        
        lh = -sum(reshape(((repmat(jls',sx,1)-A_r_cumulative).^2),numSensors,sx))';

        lh = lh * 0.5 / (w_sigma^2);                

% %THE FOLLOWING IS NOT READY FOR MULTIPLE TARGET TRACKING:
% 
%     case 1 %exponential decay
%         %TODO
%     case 2 %bearing-only
%         lh = zeros(sx,1);       
%         for k = 1:numSensors
%             diffA  = abs(jls(k) - atan2(x_n(:,2)-sensorsPos(k,2),x_n(:,1)-sensorsPos(k,1))); %abs value of the difference in the exponent of the likelihood
%             diffA(diffA > pi) = 2*pi - diffA(diffA > pi);
%             lh = lh - diffA.^2;
% %             lh = lh - (jls(k) - atan2(x_n(:,2)-sensorsPos(k,2),x_n(:,1)-sensorsPos(k,1))).^2;   %OLD AND NOT WORKING SOLUTION -> just for illustration...         
%         end
%         
%         lh = lh * 0.5 / (w_sigma^2);
end    
%<-------------------------------------------------------------------------
else
%Consensus:--------------------------------------------------------------->

%numJLFsums = size(parset.powersJLF,1);
%lh = ones(parset.numParticles,numJLFsums-1);

%NOTE: The exponent of the JLF is a multi-variate polynomial.

%Evaluate all monomials for all particle values:
%for rr = 2:numJLFsums
%    for variable = 1:2*parset.numTargets
%        lh(:,rr-1) = lh(:,rr-1) .* x_n(:,variable).^parset.powersJLF(rr,variable);
%    end
%end

%%%% optimized code %%%
lh = exp(real(log(x_n)*(parset.powersJLF(2:end,:))')); %powers(2:end) since the first row are all zeros
% lh = exp(log(x_n)*(parset.powersJLF(2:end,:))'); %powers(2:end) since the first row are all zeros
% lh = real(lh); %due to log(x) for non-positive x

%NOTE: Non-positive x_n are very rare. They occure only when target is
%close to left or bottom border of simulation area. In that case it can
%happen that some particles are sampled out of bounds. In that case the
%above code is only approximate, but it shouldn't have any signifficant 
%effects on the final performance.
%%%%%%%%%%%%%%%%%%%%%%%

%Linear combination of the monomials: (coefficients are given by jls
%computed using consensus...)
lh = lh*jls.';

lh = -lh./(2*parset.w_sigma^2);
%<-------------------------------------------------------------------------
end

lh = lh - max(lh);
lh = exp(lh);

end