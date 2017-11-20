% %THE FOLLOWING IS NOT READY FOR MULTIPLE TARGET TRACKING:
% function lh = singleLikelihood(x_n,sensor,z)
% %Computes likelihood values for each input sample for the given sensor and
% %its measurement.
% %
% %x_n - numSamples x 2 vector containing the samples where the likelihood
% %is to be evaluated
% %sensor - index of sensor for which likelihood is to be evaluated
% %z - measurement of the sensor
% %
% %lh - numSamples x 1 vector of likelihood values evaluated at x_n
% %(likelihood is evaluated up to a multiplicative constant factor)
% 
% global parset;
% 
% sensorPos = parset.sensorsPos(sensor,:);%to speed up execution
% 
% switch parset.measModel
%     case 0 %inverse decay
%         sx = size(x_n,1);%number of samples
% %         lh = zeros(sx,1);        
% %         A = parset.A;%to speed up execution
% %         
% %         invPow = parset.invPow;%to speed up execution        
% %         for s = 1:sx
% %                 r = norm(x_n(s,:) - sensorPos);
% %                 lh(s) = -(z - A/(r^invPow)).^2;
% %         end
% %         lh = lh * 0.5 / (parset.w_sigma^2);
%                 
%         r = x_n - repmat(sensorPos,sx,1);
%         r = hypot(r(:,1),r(:,2));
%         lh = -(z - parset.A./(r.^parset.invPow)).^2;
% 
%         lh = lh * 0.5 / (parset.w_sigma^2);
%         
%     case 1 %exponential decay
%         %TODO
%     case 2 %bearing-only
%         lh = -(0.5/(parset.w_sigma^2))*(z - atan2(x_n(:,2)-sensorPos(2),x_n(:,1)-sensorPos(1))).^2;    
% end    
% 
% lh = lh - max(lh);
% lh = exp(lh);
% 
% end