function jointLikelihoodSums = averageConsensus(z_n)

global parset;
global particlesApprox; %contains particle positions
%used as "sample points" in least squares approximation and to calculate
%reference points for Taylor approximation

persistent P L targetPosIndx;

consensusNumOfIter = parset.consensusNumOfIter;

G   = parset.adjacencyMat;          %adjacenecy
D   = sum(G,2);                     %degrees
N   = size(G,1);

switch parset.consensusWeightModel(1)
    case 0 %constant step (uniform weights)
        eta = parset.consensusWeightModel(2);     %step size

        if (isempty(P))
            L = sparse(diag(D))-G;                            %Laplacian matrix
            P = speye(parset.numSensors)-eta*L;         %Perron matrix
        end
    case 1 %maximum degree weights
        eta = 1/max(D);
        if (isempty(P))
            L = sparse(diag(D))-G;                            %Laplacian matrix
            P = speye(parset.numSensors)-eta*L;         %Perron matrix
        end
    case 2 %local degree weights
        if (isempty(P))
            P = zeros(N);
            for i=1:N
                for j=(i+1):N
                    if (G(i,j)==1)
                        P(i,j) = 1/(max(D(i),D(j)));
                    end
                end
            end
            P = P+P';
            dg = diag(1-sum(P,2));
            P = sparse(P+dg);
        end
    case 3 %metropolis weights
        if (isempty(P))
            P = zeros(N);
            for i=1:N
                for j=(i+1):N
                    if (G(i,j)==1)
                        P(i,j) = 1/(1+max(D(i),D(j)));
                    end
                end
            end
            P = P+P';
            dg = diag(1-sum(P,2));
            P = sparse(P+dg);
        end
end


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
%     
%         p = parset.invPow;
%         s1 = parset.sensorsPos(:,1);
%         s2 = parset.sensorsPos(:,2);
%         A = parset.A;
%         
%         z = z_n;
%         x1p = x_prime(:,1);
%         x2p = x_prime(:,2);
%         
%         h_0_0 = A./((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2);
%         h_0_1 = (A*p*(2*s2 - 2*x2p))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 1));
%         h_1_0 = (A*p*(2*s1 - 2*x1p))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 1));
%         h_0_2 = (A*p*(2*s2 - 2*x2p).^2*(p/2 + 1))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 2)) - (A*p)./((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 1);
%         h_1_1 = (A*p*(2*s1 - 2*x1p).*(2*s2 - 2*x2p).*(p/2 + 1))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 2));
%         h_2_0 = (A*p*(2*s1 - 2*x1p).^2.*(p/2 + 1))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 2)) - (A*p)./((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 1);         
%   
%         x1 = h_2_0.^2/4;
%         x2 = h_1_1.*h_2_0; 
%         x3 = -h_2_0.*(h_2_0.*x1p - h_1_0 + h_1_1.*x2p); 
%         x4 = h_1_1.^2 + (h_0_2.*h_2_0)/2; 
%         x5 = h_0_1.*h_2_0 + 2*h_1_0.*h_1_1 - 2*h_1_1.^2.*x2p - 3*h_1_1.*h_2_0.*x1p - h_0_2.*h_2_0.*x2p; 
%         x6 = h_0_0.*h_2_0 - h_2_0.*z + h_1_0.^2 + (3*h_2_0.^2.*x1p.^2)/2 + h_1_1.^2.*x2p.^2 + (h_0_2.*h_2_0.*x2p.^2)/2 - 3*h_1_0.*h_2_0.*x1p - h_0_1.*h_2_0.*x2p - 2*h_1_0.*h_1_1.*x2p + 3*h_1_1.*h_2_0.*x1p.*x2p;
%         x7 = h_0_2.*h_1_1; 
%         x8 = 2*h_0_1.*h_1_1 + h_0_2.*h_1_0 - 2*h_1_1.^2.*x1p - h_0_2.*h_2_0.*x1p - 3*h_0_2.*h_1_1.*x2p;
%         x9 = 2*h_0_0.*h_1_1 + 2*h_0_1.*h_1_0 - 2*h_1_1.*z + 3*h_1_1.*h_2_0.*x1p.^2 + 3*h_0_2.*h_1_1.*x2p.^2 + 4*h_1_1.^2.*x1p.*x2p - 2*h_0_1.*h_2_0.*x1p - 4*h_1_0.*h_1_1.*x1p - 4*h_0_1.*h_1_1.*x2p - 2*h_0_2.*h_1_0.*x2p + 2*h_0_2.*h_2_0.*x1p.*x2p;
%         x10 = -(h_2_0.*x1p - h_1_0 + h_1_1.*x2p).*(h_2_0.*x1p.^2 + 2*h_1_1.*x1p.*x2p - 2*h_1_0.*x1p + h_0_2.*x2p.^2 - 2*h_0_1.*x2p + 2*h_0_0 - 2*z);
%         x11 = h_0_2.^2/4; 
%         x12 = -h_0_2.*(h_1_1.*x1p - h_0_1 + h_0_2.*x2p); 
%         x13 = h_0_0.*h_0_2 - h_0_2.*z + h_0_1.^2 + h_1_1.^2.*x1p.^2 + (3*h_0_2.^2.*x2p.^2)/2 + (h_0_2.*h_2_0.*x1p.^2)/2 - 2*h_0_1.*h_1_1.*x1p - h_0_2.*h_1_0.*x1p - 3*h_0_1.*h_0_2.*x2p + 3*h_0_2.*h_1_1.*x1p.*x2p;
%         x14 = -(h_1_1.*x1p - h_0_1 + h_0_2.*x2p).*(h_2_0.*x1p.^2 + 2*h_1_1.*x1p.*x2p - 2*h_1_0.*x1p + h_0_2.*x2p.^2 - 2*h_0_1.*x2p + 2*h_0_0 - 2*z); 
%         
%         jointLikelihoodSums = (P^consensusNumOfIter)*[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14];        
%         jointLikelihoodSums = parset.numSensors*jointLikelihoodSums; %get sum from the average
% 
%     case 1 %exponential decay
%         %TODO
%     case 2 %bearing-only
%         %TODO
%         
% end

    case 1 %Least squares approximation
    
    %NOTE: least squares approx. of 2nd order used
    %NOTE: only inverse decay measurement model is implemented!
            
%     alpha = zeros(parset.numSensors,6); %2nd order LS approx... => each row contains 6 coeffs
%     
%     for k=1:parset.numSensors %calculate the LS coeffs (alpha) for each sensor
%         
%         particlesApprox_k = particlesApprox(:,:,k);
%         %NOTE: The line below is maybe not neccessary. Commented to speed up...
%         %particlesApprox_k = unique(particlesApprox_k,'rows'); %remove repeating particles
%         
%         %inverse decay measurement model:
%         hVector = parset.A ./ ((hypot(particlesApprox_k(:,1)-parset.sensorsPos(k,1),particlesApprox_k(:,2)-parset.sensorsPos(k,2))).^(parset.invPow));
%         hVector(~isfinite(hVector)) = parset.A;
%         
%         x1 = particlesApprox_k(:,1);
%         x2 = particlesApprox_k(:,2);
%         polyTermsMat = [ones(size(x1)) x1 x2 x1.*x2 x1.^2 x2.^2];
%         alpha(k,:) = polyTermsMat\hVector; %polynomial coefficients obtained by least squares
%     end
% 
%     z = z_n;
%    
%     tmp1 = alpha(:,5).^2;
%     tmp2 = 2*alpha(:,4).*alpha(:,5);
%     tmp3 = 2*alpha(:,2).*alpha(:,5);
%     tmp4 = alpha(:,4).^2 + 2*alpha(:,6).*alpha(:,5);
%     tmp5 = 2*(alpha(:,3).*alpha(:,5) + alpha(:,2).*alpha(:,4));
%     tmp6 = alpha(:,2).^2 + 2*alpha(:,1).*alpha(:,5) - 2*alpha(:,5).*z;
%     tmp7 = 2*alpha(:,6).*alpha(:,4);
%     tmp8 = 2*(alpha(:,3).*alpha(:,4) + alpha(:,2).*alpha(:,6));
%     tmp9 = 2*(alpha(:,1).*alpha(:,4) + alpha(:,3).*alpha(:,2) - alpha(:,4).*z);
%     tmp10 = 2*alpha(:,2).*(alpha(:,1) - z);
%     tmp11 = alpha(:,6).^2;
%     tmp12 = 2*alpha(:,3).*alpha(:,6);
%     tmp13 = alpha(:,3).^2 + 2*alpha(:,1).*alpha(:,6) - 2*alpha(:,6).*z;
%     tmp14 = 2*alpha(:,3).*(alpha(:,1) - z);
%     
%     jointLikelihoodSums = (P^consensusNumOfIter)*[tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14];        
%     jointLikelihoodSums = parset.numSensors*jointLikelihoodSums; %get sum from the average

     switch parset.measModel
         case 0 %inverse decay
            numJLFsums = size(parset.powersJLF,1);
            tmp = zeros(parset.numSensors,numJLFsums-1);
            if (isempty(targetPosIndx))
                targetPosIndx = sparse(repmat([1 1 0 0],1,parset.numTargets));
            end
          
            for k=1:parset.numSensors
                x = particlesApprox(:,targetPosIndx==1,k); %pick only the positions...

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

                tmp(k,:) = multiVarConv(alpha);

            end

            jointLikelihoodSums = (P^consensusNumOfIter)*tmp;        
            jointLikelihoodSums = parset.numSensors*jointLikelihoodSums; %get sum from the average

         case 1 %exponential decay
           %TODO
         case 2 %bearing-only
            %TODO
     end
         
end

end