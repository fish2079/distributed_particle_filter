function jointLikelihoodSums = quantizedCensi(z_n)%,x_prime)

global parset;
global particlesApprox; %contains particle positions
%used as "sample points" in least squares approximation and to calculate
%reference points for Taylor approximation

persistent tmpA tmpD I J D G targetPosIndx polyTermsMat; %optimized


consensusNumOfIter = parset.consensusNumOfIter;

switch parset.consensusWeightModel(1)
    case 0
        eta = parset.consensusWeightModel(2);              %step size
        roundType = parset.consensusAlgorithm(2);
        wordLength = parset.consensusAlgorithm(3);
        fractionLength = parset.consensusAlgorithm(4);
        formatType = parset.consensusAlgorithm(5);
        N   = parset.numSensors;
        if (isempty(G))
            G   = parset.adjacencyMat;          %adjacenecy
            D   = sparse(diag(sum(G,2)));               %degree matrix
        end
        Delta = max(max(D));
    otherwise
        error('Censi: Weight model must be set to Constant in parset.txt!');
end
        

%L = diag(D)-A;                      %Laplacian matrix
%P = eye(parset.numSensors)-eta*L;   %Perron matrix

switch parset.polynomialApprox
    case 0 %Taylor approximation
    switch parset.measModel
        case 0 %inverse decay
           %%
    %         p = parset.invPow;
    %         s1 = parset.sensorsPos(:,1);
    %         s2 = parset.sensorsPos(:,2);
    %         A = parset.A;
    %                 
    %         z = z_n;
    %         x1p = x_prime(:,1);
    %         x2p = x_prime(:,2);
    %         
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
    %         x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12;x13;x14];
    %         
    %         c = zeros(14*N,1);
    %         
    %         tmpD = eye(14*N)-eta/Delta*kron(eye(14),D);
    %         tmpA = eta/Delta*kron(eye(14),G);
    %         
    %         for i=1:consensusNumOfIter
    %             y = quant(x-c,0);
    %             c = c+(y-x);
    %             x = tmpD*x+tmpA*y;
    %         end
    %         
    %         jointLikelihoodSums = zeros(N,14);
    %         jointLikelihoodSums(1:end) = N*x;
    %         
            %%%%%%%%%
% 
%             z = z_n;        
% 
%             s1 = parset.sensorsPos(:,1);
%             s2 = parset.sensorsPos(:,2);
% 
%             x1p = x_prime(:,1);
%             x2p = x_prime(:,2);
% 
%             h_0_0 = atan2((x2p - s2),(x1p - s1));
%             h_0_1 = -1./((s1 - x1p).*((s2 - x2p).^2./(s1 - x1p).^2 + 1));
%             h_1_0 = (s2 - x2p)./((s1 - x1p).^2.*((s2 - x2p).^2./(s1 - x1p).^2 + 1));
%             h_0_2 = -(2*s2 - 2*x2p)./((s1 - x1p).^3.*((s2 - x2p).^2./(s1 - x1p).^2 + 1).^2);
%             h_1_1 = ((2*s2 - 2*x2p).*(s2 - x2p))./((s1 - x1p).^4.*((s2 - x2p).^2./(s1 - x1p).^2 + 1).^2) - 1./((s1 - x1p).^2.*((s2 - x2p).^2./(s1 - x1p).^2 + 1));
%             h_2_0 = (2*(s2 - x2p))./((s1 - x1p).^3.*((s2 - x2p).^2./(s1 - x1p).^2 + 1)) - (2*(s2 - x2p).^3)./((s1 - x1p).^5.*((s2 - x2p).^2./(s1 - x1p).^2 + 1).^2);            
% 
%             x1 = h_2_0.^2/4;
%             x2 = h_1_1.*h_2_0;
%             x3 = -h_2_0.*(h_2_0.*x1p - h_1_0 + h_1_1.*x2p);
%             x4 = h_1_1.^2 + (h_0_2.*h_2_0)/2;
%             x5 = h_0_1.*h_2_0 + 2*h_1_0.*h_1_1 - 2*h_1_1.^2.*x2p - 3*h_1_1.*h_2_0.*x1p - h_0_2.*h_2_0.*x2p;
%             x6 = h_0_0.*h_2_0 - h_2_0.*z + h_1_0.^2 + (3*h_2_0.^2.*x1p.^2)/2 + h_1_1.^2.*x2p.^2 + (h_0_2.*h_2_0.*x2p.^2)/2 - 3*h_1_0.*h_2_0.*x1p - h_0_1.*h_2_0.*x2p - 2*h_1_0.*h_1_1.*x2p + 3*h_1_1.*h_2_0.*x1p.*x2p;
%             x7 = h_0_2.*h_1_1;
%             x8 = 2*h_0_1.*h_1_1 + h_0_2.*h_1_0 - 2*h_1_1.^2.*x1p - h_0_2.*h_2_0.*x1p - 3*h_0_2.*h_1_1.*x2p;
%             x9 = 2*h_0_0.*h_1_1 + 2*h_0_1.*h_1_0 - 2*h_1_1.*z + 3*h_1_1.*h_2_0.*x1p.^2 + 3*h_0_2.*h_1_1.*x2p.^2 + 4*h_1_1.^2.*x1p.*x2p - 2*h_0_1.*h_2_0.*x1p - 4*h_1_0.*h_1_1.*x1p - 4*h_0_1.*h_1_1.*x2p - 2*h_0_2.*h_1_0.*x2p + 2*h_0_2.*h_2_0.*x1p.*x2p;
%             x10 = -(h_2_0.*x1p - h_1_0 + h_1_1.*x2p).*(h_2_0.*x1p.^2 + 2*h_1_1.*x1p.*x2p - 2*h_1_0.*x1p + h_0_2.*x2p.^2 - 2*h_0_1.*x2p + 2*h_0_0 - 2*z);
%             x11 = h_0_2.^2/4;
%             x12 = -h_0_2.*(h_1_1.*x1p - h_0_1 + h_0_2.*x2p);
%             x13 = h_0_0.*h_0_2 - h_0_2.*z + h_0_1.^2 + h_1_1.^2.*x1p.^2 + (3*h_0_2.^2.*x2p.^2)/2 + (h_0_2.*h_2_0.*x1p.^2)/2 - 2*h_0_1.*h_1_1.*x1p - h_0_2.*h_1_0.*x1p - 3*h_0_1.*h_0_2.*x2p + 3*h_0_2.*h_1_1.*x1p.*x2p;
%             x14 = -(h_1_1.*x1p - h_0_1 + h_0_2.*x2p).*(h_2_0.*x1p.^2 + 2*h_1_1.*x1p.*x2p - 2*h_1_0.*x1p + h_0_2.*x2p.^2 - 2*h_0_1.*x2p + 2*h_0_0 - 2*z); 
% 
%             x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12;x13;x14];
% 
%             c = zeros(14*N,1);
% 
%             tmpD = eye(14*N)-eta/Delta*kron(eye(14),D);
%             tmpA = eta/Delta*kron(eye(14),G);
% 
%             for i=1:consensusNumOfIter
%                 y = quant(x-c,0);
%                 c = c+(y-x);
%                 x = tmpD*x+tmpA*y;
%             end
% 
%             jointLikelihoodSums = zeros(N,14);
%             jointLikelihoodSums(1:end) = N*x;
        case 1  %exponential decay
           %TODO
        case 2 %bearing-only
            %TODO
    end
    
    case 1 %Least squares approximation
        switch parset.measModel
            case 0 %inverse decay
                numJLFsums = size(parset.powersJLF,1);
                tmp = zeros(N,numJLFsums-1);
                if (isempty(targetPosIndx))
                    targetPosIndx = sparse(repmat([1 1 0 0],1,parset.numTargets));
                end
                
                for k=1:parset.numSensors
                    x = particlesApprox(:,targetPosIndx==1,k); %pick only the positions...
                
                    %Inverse decay measurement model:
                    hVector = zeros(parset.numParticles,1);
                    for tt = 0:parset.numTargets-1
                        A_r = parset.A ./ ((hypot(x(:,tt*2+1)-parset.sensorsPos(k,1),x(:,tt*2+2)-parset.sensorsPos(k,2))).^(parset.invPow));
                        A_r(A_r > parset.A) = parset.A;
                        hVector = hVector + A_r;
                    end
                    
%                     polyTermsMat = ones(parset.numParticles,size(parset.powersMeasFunc,1)); %matrix of monomials evaluated at particles
% 
%                     for monom = 2:size(parset.powersMeasFunc,1) %the first monomial is constant equal to 1 so we skip it...
%                         for variable = 1:2*parset.numTargets
%                             polyTermsMat(:,monom) = polyTermsMat(:,monom) .* x(:,variable).^parset.powersMeasFunc(monom,variable); 
%                         end
%                     end 

                    %%%% optimized version %%%% only faster than for loop(above) if
                    %%%% parset.powersMeasFunc is sparse
                    polyTermsMat = exp(log(x)*(parset.powersMeasFunc)');
                    polyTermsMat = real(polyTermsMat); %due to log(x) for negative x
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%

                    alpha = polyTermsMat\hVector; %polynomial coefficients obtained by least squares

                    alpha(1) = alpha(1) - z_n(k);
                    tmp(k,:) = multiVarConv(alpha);

                end
                
                M = N*(numJLFsums-1);
                c = zeros(M,1);
                x = tmp(:);
           
                if (isempty(tmpA))
                    I = speye(M);
                    J = speye(numJLFsums-1);
                
                    tmpD = I-eta/Delta*kron(J,D);
                    tmpA = eta/Delta*kron(J,G);
                end
                
                for i=1:consensusNumOfIter
                    y = quant(x-c,roundType,wordLength,fractionLength,formatType);
                    c = c+(y-x);
                    x = tmpD*x+tmpA*y;
                end
                
                jointLikelihoodSums = zeros(N,(numJLFsums-1));
                jointLikelihoodSums(1:end) = N*x;

%                L = D-G;
%                P = eye(parset.numSensors)-eta/Delta*L;
%                jointLikelihoodSums2 = (P^consensusNumOfIter)*tmp;
%                jointLikelihoodSums2 = parset.numSensors*jointLikelihoodSums2; %get sum from the average

            case 1  %exponential decay
                %TODO
            case 2 %bearing-only
                %TODO
        end
end

end

function Y = quant(X,roundType,wordLength,fractionLength,formatType)

persistent  dispRange;
qAux = 0;
switch roundType
    case 0
        if (formatType == 0)
            q = quantizer('mode','fixed','roundmode','round','overflowmode','saturate','format',[wordLength fractionLength]);
            %range: -(2^(wordLength-1))/2^fractionLength -- (2^(wordLength-1)-1)/2^fractionLength
            Y = quantize(q,X);
            qAux = 1;            
        end
        if(formatType == 1)
            q = quantizer('mode','float','roundmode','round','format',[wordLength fractionLength]);
            Y = quantize(q,X);
            qAux = 1;
        end
        if(formatType == 2)     %ex.: X=10.439, wordLength=8, fractionLength=2 => Y=8.44
            if (wordLength<Inf)
                X = X-fix(X/(10^(wordLength)))*10^(wordLength);
            end
            if (fractionLength<Inf)
                Y = (round(X*10^fractionLength))/10^fractionLength;
            else
                Y=X;
            end
        end
        
        %Y = round(X);
    case 1
        if (formatType == 0)
            q = quantizer('mode','fixed','roundmode','ceil','overflowmode','saturate','format',[wordLength fractionLength]);
            Y = quantize(q,X);
            qAux = 1;
            %range: -(2^(wordLength-1))/2^fractionLength -- (2^(wordLength-1)-1)/2^fractionLength
        end
        if (formatType == 1)
            q = quantizer('mode','float','roundmode','ceil','format',[wordLength fractionLength]);
            Y = quantize(q,X);
            qAux = 1;
        end
        if(formatType == 2)   %ex.: X=10.439, wordLength=8, fractionLength=2 => Y=8.44
            if (wordLength<Inf)
                X = X-fix(X/(10^(wordLength)))*10^(wordLength);
            end
            if (fractionLength<Inf)
                Y = (ceil(X*10^fractionLength))/10^fractionLength;
            else
                Y = X;
            end
        end

        %Y = ceil(X);
    case 2
        if (formatType == 0)
            q = quantizer('mode','fixed','roundmode','floor','overflowmode','saturate','format',[wordLength fractionLength]);
            Y = quantize(q,X);
            qAux = 1;
            %range: -(2^(wordLength-1))/2^fractionLength -- (2^(wordLength-1)-1)/2^fractionLength
        end
        if (formatType == 1)
           q = quantizer('mode','float','roundmode','floor','format',[wordLength fractionLength]);
           Y = quantize(q,X);
           qAux = 1;
        end
        if (formatType == 2) %ex.: X=10.439, wordLength=8, fractionLength=2 => Y=8.43
           if (wordLength<Inf)
               X = X-fix(X/(10^(wordLength)))*10^(wordLength);
           end
           if (fractionLength<Inf)
               Y = (floor(X*10^fractionLength))/10^fractionLength;
           else
               Y = X;
           end
        end

        %Y = floor(X);
end
         
if (isempty(dispRange) && qAux == 1)
    dispRange = range(q);
    disp(['Quantizer range: ',num2str(dispRange)]);
end

end
