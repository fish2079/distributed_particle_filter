function jointLikelihoodSums = likelihoodConsensus_perfect(z_n)
%Function calculates the sums defining the joint likelihood for each
%sensor. The sums are computed exactly and therefore are equal for all
%sensors.
%
%z_n - vector containing measurements of all sensors
%
%jointLikelihoodSums - matrix storing in rows the consensus sums defining 
%the joint likelihood for each sensor; since we consider perfect consensus
%here, the rows of this matrix are identical (we however use entire matrix
%due to compatibility with particle filtering functions)

global parset;

global particlesApprox; %contains predicted particles
%These particles are used as "sample points" in least squares approximation
%or to calculate reference points x_prime in Taylor approximation.

switch parset.polynomialApprox
    case 0 %Taylor approximation

% %NOTE: The code for "Taylor approximation" does not support tracking of
% %multiple targets and arbitrary degree of the approximating polynomial.
% %Only one target and 2nd degree approximation is supported.
% 
% %Calculate x_prime:----------------->
% x_prime = zeros(parset.numSensors,2);
% for k=1:parset.numSensors
%     particles_positions = particlesApprox(:,1:2,k);
%     particles_weights = particlesApprox(:,5,k);
%     particles_weights = particles_weights .* singleLikelihood(particles_positions,k,z_n(k),parset);
%     particles_weights = particles_weights/sum(particles_weights);
% 
%     %Calculate target position estimate (based on local measurement only):
%     x_prime(k,:) = sum(repmat(particles_weights,1,2).*particles_positions);
% end
% %<-----------------------------------
% 
% switch parset.measModel
%     case 0 %inverse decay
%         p = parset.invPow;
%         jointLikelihoodSums = zeros(1,14);
%         A = parset.A;
%         
%         for k=1:parset.numSensors
%             
%             s1 = parset.sensorsPos(k,1);
%             s2 = parset.sensorsPos(k,2);
%             x1p = x_prime(k,1);
%             x2p = x_prime(k,2);
%             
%             z = z_n(k);
%             
%             %NOTE: Taylor approx. of 2nd order used
%             
%             h_0_0 = A/((s1 - x1p)^2 + (s2 - x2p)^2)^(p/2);
%             h_0_1 = (A*p*(2*s2 - 2*x2p))/(2*((s1 - x1p)^2 + (s2 - x2p)^2)^(p/2 + 1));
%             h_1_0 = (A*p*(2*s1 - 2*x1p))/(2*((s1 - x1p)^2 + (s2 - x2p)^2)^(p/2 + 1));
%             h_0_2 = (A*p*(2*s2 - 2*x2p)^2*(p/2 + 1))/(2*((s1 - x1p)^2 + (s2 - x2p)^2)^(p/2 + 2)) - (A*p)/((s1 - x1p)^2 + (s2 - x2p)^2)^(p/2 + 1);
%             h_1_1 = (A*p*(2*s1 - 2*x1p)*(2*s2 - 2*x2p)*(p/2 + 1))/(2*((s1 - x1p)^2 + (s2 - x2p)^2)^(p/2 + 2));
%             h_2_0 = (A*p*(2*s1 - 2*x1p)^2*(p/2 + 1))/(2*((s1 - x1p)^2 + (s2 - x2p)^2)^(p/2 + 2)) - (A*p)/((s1 - x1p)^2 + (s2 - x2p)^2)^(p/2 + 1);         
%  
%             jointLikelihoodSums(1) = jointLikelihoodSums(1) + h_2_0^2/4;
%             jointLikelihoodSums(2) = jointLikelihoodSums(2) + h_1_1*h_2_0;
%             jointLikelihoodSums(3) = jointLikelihoodSums(3) + -h_2_0*(h_2_0*x1p - h_1_0 + h_1_1*x2p);
%             jointLikelihoodSums(4) = jointLikelihoodSums(4) + h_1_1^2 + (h_0_2*h_2_0)/2;
%             jointLikelihoodSums(5) = jointLikelihoodSums(5) + h_0_1*h_2_0 + 2*h_1_0*h_1_1 - 2*h_1_1^2*x2p - 3*h_1_1*h_2_0*x1p - h_0_2*h_2_0*x2p;
%             jointLikelihoodSums(6) = jointLikelihoodSums(6) + h_0_0*h_2_0 - h_2_0*z + h_1_0^2 + (3*h_2_0^2*x1p^2)/2 + h_1_1^2*x2p^2 + (h_0_2*h_2_0*x2p^2)/2 - 3*h_1_0*h_2_0*x1p - h_0_1*h_2_0*x2p - 2*h_1_0*h_1_1*x2p + 3*h_1_1*h_2_0*x1p*x2p;
%             jointLikelihoodSums(7) = jointLikelihoodSums(7) + h_0_2*h_1_1;
%             jointLikelihoodSums(8) = jointLikelihoodSums(8) + 2*h_0_1*h_1_1 + h_0_2*h_1_0 - 2*h_1_1^2*x1p - h_0_2*h_2_0*x1p - 3*h_0_2*h_1_1*x2p;
%             jointLikelihoodSums(9) = jointLikelihoodSums(9) + 2*h_0_0*h_1_1 + 2*h_0_1*h_1_0 - 2*h_1_1*z + 3*h_1_1*h_2_0*x1p^2 + 3*h_0_2*h_1_1*x2p^2 + 4*h_1_1^2*x1p*x2p - 2*h_0_1*h_2_0*x1p - 4*h_1_0*h_1_1*x1p - 4*h_0_1*h_1_1*x2p - 2*h_0_2*h_1_0*x2p + 2*h_0_2*h_2_0*x1p*x2p;
%             jointLikelihoodSums(10) = jointLikelihoodSums(10) + -(h_2_0*x1p - h_1_0 + h_1_1*x2p)*(h_2_0*x1p^2 + 2*h_1_1*x1p*x2p - 2*h_1_0*x1p + h_0_2*x2p^2 - 2*h_0_1*x2p + 2*h_0_0 - 2*z);
%             jointLikelihoodSums(11) = jointLikelihoodSums(11) + h_0_2^2/4;
%             jointLikelihoodSums(12) = jointLikelihoodSums(12) + -h_0_2*(h_1_1*x1p - h_0_1 + h_0_2*x2p);
%             jointLikelihoodSums(13) = jointLikelihoodSums(13) + h_0_0*h_0_2 - h_0_2*z + h_0_1^2 + h_1_1^2*x1p^2 + (3*h_0_2^2*x2p^2)/2 + (h_0_2*h_2_0*x1p^2)/2 - 2*h_0_1*h_1_1*x1p - h_0_2*h_1_0*x1p - 3*h_0_1*h_0_2*x2p + 3*h_0_2*h_1_1*x1p*x2p;
%             jointLikelihoodSums(14) = jointLikelihoodSums(14) + -(h_1_1*x1p - h_0_1 + h_0_2*x2p)*(h_2_0*x1p^2 + 2*h_1_1*x1p*x2p - 2*h_1_0*x1p + h_0_2*x2p^2 - 2*h_0_1*x2p + 2*h_0_0 - 2*z); 
%             
%         end
% 
%     case 1 %exponential decay
%         %TODO
%     case 2 %bearing-only
%         %TODO
%         
% end

    case 1 %Least squares approximation

      switch parset.measModel
         case 0 %inverse decay
            numJLFsums = size(parset.powersJLF,1);

            jointLikelihoodSums = zeros(1,numJLFsums-1);

            for k=1:parset.numSensors

                x = particlesApprox(:,repmat([1 1 0 0],1,parset.numTargets)==1,k); %pick only the positions...

                %Inverse decay measurement model:
                hVector = zeros(parset.numParticles,1);
                if parset.trackLossDetOnOff > 1
                    A_r_detection = zeros(parset.numParticles,parset.numTargets);
                end
                for tt = 0:parset.numTargets-1
                    A_r = parset.A ./ ((hypot(x(:,tt*2+1)-parset.sensorsPos(k,1),x(:,tt*2+2)-parset.sensorsPos(k,2))).^(parset.invPow));
                    A_r(A_r > parset.A/(parset.d0^(parset.invPow))) = parset.A/(parset.d0^(parset.invPow));
                    hVector = hVector + A_r;
                    if parset.trackLossDetOnOff > 1
                        A_r_detection(:,tt+1) = A_r;
                    end
                end

%                 polyTermsMat = ones(parset.numParticles,size(parset.powersMeasFunc,1)); %matrix of monomials evaluated at particles
%                 
%                 for monom = 2:size(parset.powersMeasFunc,1) %the first monomial is constant equal to 1 so we skip it...
%                     for variable = 1:2*parset.numTargets
%                         polyTermsMat(:,monom) = polyTermsMat(:,monom) .* x(:,variable).^parset.powersMeasFunc(monom,variable); 
%                     end
%                 end  

%                 %%%% optimized version %%%% only faster than for loop(above) if
%                 %%%% parset.powersMeasFunc is sparse
%                 polyTermsMat = exp(real(log(x)*(parset.powersMeasFunc)'));
%                 %polyTermsMat = real(polyTermsMat);  %due to log(x) for negative x
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                 alpha = polyTermsMat\hVector; %polynomial coefficients obtained by least squares
%                 %alpha = pinv(polyTermsMat)*hVector; %polynomial coefficients obtained by least squares

%                 polymodel = polyfitn(x,hVector,parset.powersMeasFunc);
%                 alpha = polymodel.Coefficients;

%                 polyTermsMat = ones(parset.numParticles,size(parset.powersMeasFunc,1)); %matrix of monomials evaluated at particles
%                 for tt = 0:parset.numTargets-1
%                     polyTermsMat = polyTermsMat.*(bsxfun(@power,x(:,2*tt+1),parset.powersMeasFunc(:,2*tt+1)').*bsxfun(@power,x(:,2*tt+2),parset.powersMeasFunc(:,2*tt+2)'));
%                 end
%                 alpha=polyTermsMat\hVector; %polynomial coefficients obtained by least squares

                if parset.trackLossDetOnOff == 2 && sum(sum(A_r_detection > (1/(parset.trackLossDetTrshld^parset.invPow))*parset.A))>0
                    %diregard sensors that are too close to one of the
                    %targets:
                    alpha = zeros(size(parset.powersMeasFunc,1),1);
%                     disp('ZEROS');
                    %NOTE: Setting alphas to 0 means that we approximate
                    %the measurement function of this sensor by a 0
                    %polynomial. This is the measurement function of a
                    %sensor that is in an infinite distance from the
                    %targets. => we disregard this sensor => this improves
                    %numerical stability of particle filterig
                else
%                     polyTermsMat = ones(parset.numParticles,size(parset.powersMeasFunc,1)); %matrix of monomials evaluated at particles
%                     for tt = 0:parset.numTargets-1
%                         polyTermsMat = polyTermsMat.*(bsxfun(@power,x(:,2*tt+1),parset.powersMeasFunc(:,2*tt+1)').*bsxfun(@power,x(:,2*tt+2),parset.powersMeasFunc(:,2*tt+2)'));
%                     end

                    polyTermsMat = exp(real(log(x)*(parset.powersMeasFunc)'));
                    
                    [Q,R]=qr(polyTermsMat,0);
                    alpha=R\(Q'*hVector); %polynomial coefficients obtained by least squares
                    
%                     alpha=polyTermsMat\hVector; %polynomial coefficients obtained by least squares

                    alpha(1) = alpha(1) - z_n(k);
                end

                %NOTE: Solving Ax=b: If A has more rows than columns and is not 
                %of full rank, then the overdetermined least squares problem 
                %"min norm(A*x-b)" does not have a unique solution. Two of the 
                %infinitely many solutions are x = pinv(A)*b and y = A\b. These two
                %are distinguished by the facts that norm(x) is smaller than the 
                %norm of any other solution and that y has the fewest possible 
                %nonzero components. Both of these are exact solutions in the sense
                %that norm(A*x-b) and norm(A*y-b) are on the order of roundoff 
                %error.  (Check the docs of pinv for details.)

                jointLikelihoodSums = jointLikelihoodSums + multiVarConv(alpha); %summation over all sensors

            end
        case 1 %exponential decay
        %TODO
        case 2 %bearing-only
        %TODO
      end
end

jointLikelihoodSums = repmat(jointLikelihoodSums,parset.numSensors,1); %contains likelihood sums for each sensor (stored in rows)
%NOTE: Due to "perfect" consensus, the likelihood sums are equal for all
%sensors. In realistic consensus this is not the case.

end