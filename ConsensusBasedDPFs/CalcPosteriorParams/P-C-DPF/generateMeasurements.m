function z_n = generateMeasurements(x_n)
%Determine measurements for all sensors and current state x_n.
%
%x_n - current state (contains position of each target; velocities are not included in this vector)
%
%z_n - numSensors x 1 vector containing measurements for each sensor

global parset;

z_n = zeros(parset.numSensors,1);

for k=1:parset.numSensors
        
    switch parset.measModel
        case 0 %inverse decay
            for tt = 0:parset.numTargets-1 %the measurement is a sum of "received powers" from all targets... (plus noise)
                r = norm(x_n(tt*2+1:(tt+1)*2) - parset.sensorsPos(k,:)'); %distance of the sensor from the target
                if r >=parset.d0 %if the distance is greater than d0, we calculate the "received power" as usual... (as written in our paper)
                    z_n(k) = z_n(k) + parset.A / (r^parset.invPow);
                else %if the distance is smaller than d0...
                    z_n(k) = z_n(k) + parset.A / (parset.d0^parset.invPow); %set the "received power" to the value of the source amplitute measured at the distance d0
                end
            end
            z_n(k) = z_n(k) + parset.w_sigma*randn(1); %add noise...

%         %THE FOLLOWING IS NOT READY FOR MULTIPLE TARGET TRACKING:
%         case 1 %exponential decay
%             r = norm(x_n - parset.sensorsPos(k,:)'); %distance of the sensor from the target
%             z_n(k) = parset.A*exp(-parset.lambda*r) + parset.w_sigma*randn(1);
%         case 2 %bearing-only
%             z_n(k) = atan2(x_n(2) - parset.sensorsPos(k,2),x_n(1) - parset.sensorsPos(k,1)) + parset.w_sigma*randn(1);
%             if z_n(k) > pi
%                 z_n(k) = -pi + mod(z_n(k),pi);
%             end
%             if z_n(k) < -pi
%                 z_n(k) = pi - mod(-z_n(k),pi);
%             end
    end
end

end