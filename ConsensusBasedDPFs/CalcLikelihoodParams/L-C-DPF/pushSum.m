function jointLikelihoodSums = pushSum(z_n,x_prime)

global parset;

consensusNumOfIter = parset.consensusNumOfIter;

N   = parset.numSensors;
G   = parset.adjacencyMat;          

if verLessThan('matlab','7.9')
    ver = 0;
else
    ver = 1;
end

switch parset.measModel
    case 0
%         p = parset.invPow;
%         %jointLikelihoodSums = zeros(N,8);
%         s1 = parset.sensorsPos(:,1);
%         s2 = parset.sensorsPos(:,2);
% 
%         Q1 = z_n-(parset.A+(parset.A*p)/2)*q_prime.^(-p/2);
%         Q2 = parset.A * p/2 * q_prime.^(-p/2 - 1);
%         c1 = -2*Q1.*Q2 ;
%         c2 = Q2.^2;
%         
%         x1 = 2*c1.*s1 + c2.*(-4*s1.^3 - 4*s1.*(s2.^2));
%         x2 = 2*c1.*s2 + c2.*(-4*s2.^3 - 4*s2.*(s1.^2));
%         x3 = -c1 + c2.*(6*s1.^2 + 2*s2.^2);
%         x4 = -c1 + c2.*(6*s2.^2 + 2*s1.^2);
%         x5 = 8*c2.*s1.*s2;
%         x6 = - 4*c2.*s1;
%         x7 = - 4*c2.*s2;
%         x8 = c2;  
%         
%         
%         x = [x1,x2,x3,x4,x5,x6,x7,x8];
%         s = x;
%         w = [ones(1,8);zeros(N-1,8)]; %ones(N,8); if we want to compute average
%         %P = eye(N);
%         for i=1:timeStop
%             tmpA = selectNeighbor(A,0,ver);
%             %P = 0.5*P*tmpA';
%             s = 0.5*tmpA'*s;
%             w = 0.5*tmpA'*w;
%         end
%         x = s./w;
%         %jointLikelihoodSums(1:end) = x;
%         jointLikelihoodSums = reshape(x,N,8);
%        
%         
        
        %%%%%%%%%%%%
        
        z = z_n;        
        
        p = parset.invPow;
        s1 = parset.sensorsPos(:,1);
        s2 = parset.sensorsPos(:,2);
        A = parset.A;
        
        x1p = x_prime(:,1);
        x2p = x_prime(:,2);
        
        h_0_0 = A./((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2);
        h_0_1 = (A*p*(2*s2 - 2*x2p))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 1));
        h_1_0 = (A*p*(2*s1 - 2*x1p))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 1));
        h_0_2 = (A*p*(2*s2 - 2*x2p).^2*(p/2 + 1))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 2)) - (A*p)./((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 1);
        h_1_1 = (A*p*(2*s1 - 2*x1p).*(2*s2 - 2*x2p).*(p/2 + 1))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 2));
        h_2_0 = (A*p*(2*s1 - 2*x1p).^2.*(p/2 + 1))./(2*((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 2)) - (A*p)./((s1 - x1p).^2 + (s2 - x2p).^2).^(p/2 + 1);         
  
        x1 = h_2_0.^2/4;
        x2 = h_1_1.*h_2_0;
        x3 = -h_2_0.*(h_2_0.*x1p - h_1_0 + h_1_1.*x2p);
        x4 = h_1_1.^2 + (h_0_2.*h_2_0)/2;
        x5 = h_0_1.*h_2_0 + 2*h_1_0.*h_1_1 - 2*h_1_1.^2.*x2p - 3*h_1_1.*h_2_0.*x1p - h_0_2.*h_2_0.*x2p;
        x6 = h_0_0.*h_2_0 - h_2_0.*z + h_1_0.^2 + (3*h_2_0.^2.*x1p.^2)/2 + h_1_1.^2.*x2p.^2 + (h_0_2.*h_2_0.*x2p.^2)/2 - 3*h_1_0.*h_2_0.*x1p - h_0_1.*h_2_0.*x2p - 2*h_1_0.*h_1_1.*x2p + 3*h_1_1.*h_2_0.*x1p.*x2p;
        x7 = h_0_2.*h_1_1;
        x8 = 2*h_0_1.*h_1_1 + h_0_2.*h_1_0 - 2*h_1_1.^2.*x1p - h_0_2.*h_2_0.*x1p - 3*h_0_2.*h_1_1.*x2p;
        x9 = 2*h_0_0.*h_1_1 + 2*h_0_1.*h_1_0 - 2*h_1_1.*z + 3*h_1_1.*h_2_0.*x1p.^2 + 3*h_0_2.*h_1_1.*x2p.^2 + 4*h_1_1.^2.*x1p.*x2p - 2*h_0_1.*h_2_0.*x1p - 4*h_1_0.*h_1_1.*x1p - 4*h_0_1.*h_1_1.*x2p - 2*h_0_2.*h_1_0.*x2p + 2*h_0_2.*h_2_0.*x1p.*x2p;
        x10 = -(h_2_0.*x1p - h_1_0 + h_1_1.*x2p).*(h_2_0.*x1p.^2 + 2*h_1_1.*x1p.*x2p - 2*h_1_0.*x1p + h_0_2.*x2p.^2 - 2*h_0_1.*x2p + 2*h_0_0 - 2*z);
        x11 = h_0_2.^2/4;
        x12 = -h_0_2.*(h_1_1.*x1p - h_0_1 + h_0_2.*x2p);
        x13 = h_0_0.*h_0_2 - h_0_2.*z + h_0_1.^2 + h_1_1.^2.*x1p.^2 + (3*h_0_2.^2.*x2p.^2)/2 + (h_0_2.*h_2_0.*x1p.^2)/2 - 2*h_0_1.*h_1_1.*x1p - h_0_2.*h_1_0.*x1p - 3*h_0_1.*h_0_2.*x2p + 3*h_0_2.*h_1_1.*x1p.*x2p;
        x14 = -(h_1_1.*x1p - h_0_1 + h_0_2.*x2p).*(h_2_0.*x1p.^2 + 2*h_1_1.*x1p.*x2p - 2*h_1_0.*x1p + h_0_2.*x2p.^2 - 2*h_0_1.*x2p + 2*h_0_0 - 2*z); 
        
        x = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14];
        
        s = x;
        w = [ones(1,14);zeros(N-1,14)]; %ones(N,14); if we want to compute average
        
        for i=1:consensusNumOfIter
            tmpA = selectNeighbor(G,0,ver);
            %P = 0.5*P*tmpA';
            s = 0.5*tmpA'*s;
            w = 0.5*tmpA'*w;
        end
        jointLikelihoodSums = s./w;
        %jointLikelihoodSums(1:end) = x;
        %jointLikelihoodSums = reshape(x,N,14);
        
    case 1
    case 2      
        z = z_n;        
        
        s1 = parset.sensorsPos(:,1);
        s2 = parset.sensorsPos(:,2);
        
        
        x1p = x_prime(:,1);
        x2p = x_prime(:,2);
        
        h_0_0 = atan2((x2p - s2),(x1p - s1));
        h_0_1 = -1./((s1 - x1p).*((s2 - x2p).^2./(s1 - x1p).^2 + 1));
        h_1_0 = (s2 - x2p)./((s1 - x1p).^2.*((s2 - x2p).^2./(s1 - x1p).^2 + 1));
        h_0_2 = -(2*s2 - 2*x2p)./((s1 - x1p).^3.*((s2 - x2p).^2./(s1 - x1p).^2 + 1).^2);
        h_1_1 = ((2*s2 - 2*x2p).*(s2 - x2p))./((s1 - x1p).^4.*((s2 - x2p).^2./(s1 - x1p).^2 + 1).^2) - 1./((s1 - x1p).^2.*((s2 - x2p).^2./(s1 - x1p).^2 + 1));
        h_2_0 = (2*(s2 - x2p))./((s1 - x1p).^3.*((s2 - x2p).^2./(s1 - x1p).^2 + 1)) - (2*(s2 - x2p).^3)./((s1 - x1p).^5.*((s2 - x2p).^2./(s1 - x1p).^2 + 1).^2);            
        
        x1 = h_2_0.^2/4;
        x2 = h_1_1.*h_2_0;
        x3 = -h_2_0.*(h_2_0.*x1p - h_1_0 + h_1_1.*x2p);
        x4 = h_1_1.^2 + (h_0_2.*h_2_0)/2;
        x5 = h_0_1.*h_2_0 + 2*h_1_0.*h_1_1 - 2*h_1_1.^2.*x2p - 3*h_1_1.*h_2_0.*x1p - h_0_2.*h_2_0.*x2p;
        x6 = h_0_0.*h_2_0 - h_2_0.*z + h_1_0.^2 + (3*h_2_0.^2.*x1p.^2)/2 + h_1_1.^2.*x2p.^2 + (h_0_2.*h_2_0.*x2p.^2)/2 - 3*h_1_0.*h_2_0.*x1p - h_0_1.*h_2_0.*x2p - 2*h_1_0.*h_1_1.*x2p + 3*h_1_1.*h_2_0.*x1p.*x2p;
        x7 = h_0_2.*h_1_1;
        x8 = 2*h_0_1.*h_1_1 + h_0_2.*h_1_0 - 2*h_1_1.^2.*x1p - h_0_2.*h_2_0.*x1p - 3*h_0_2.*h_1_1.*x2p;
        x9 = 2*h_0_0.*h_1_1 + 2*h_0_1.*h_1_0 - 2*h_1_1.*z + 3*h_1_1.*h_2_0.*x1p.^2 + 3*h_0_2.*h_1_1.*x2p.^2 + 4*h_1_1.^2.*x1p.*x2p - 2*h_0_1.*h_2_0.*x1p - 4*h_1_0.*h_1_1.*x1p - 4*h_0_1.*h_1_1.*x2p - 2*h_0_2.*h_1_0.*x2p + 2*h_0_2.*h_2_0.*x1p.*x2p;
        x10 = -(h_2_0.*x1p - h_1_0 + h_1_1.*x2p).*(h_2_0.*x1p.^2 + 2*h_1_1.*x1p.*x2p - 2*h_1_0.*x1p + h_0_2.*x2p.^2 - 2*h_0_1.*x2p + 2*h_0_0 - 2*z);
        x11 = h_0_2.^2/4;
        x12 = -h_0_2.*(h_1_1.*x1p - h_0_1 + h_0_2.*x2p);
        x13 = h_0_0.*h_0_2 - h_0_2.*z + h_0_1.^2 + h_1_1.^2.*x1p.^2 + (3*h_0_2.^2.*x2p.^2)/2 + (h_0_2.*h_2_0.*x1p.^2)/2 - 2*h_0_1.*h_1_1.*x1p - h_0_2.*h_1_0.*x1p - 3*h_0_1.*h_0_2.*x2p + 3*h_0_2.*h_1_1.*x1p.*x2p;
        x14 = -(h_1_1.*x1p - h_0_1 + h_0_2.*x2p).*(h_2_0.*x1p.^2 + 2*h_1_1.*x1p.*x2p - 2*h_1_0.*x1p + h_0_2.*x2p.^2 - 2*h_0_1.*x2p + 2*h_0_0 - 2*z); 
        
        x = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14];
        
        s = x;
        w = [ones(1,14);zeros(N-1,14)]; %ones(N,14); if we want to compute average
        
        for i=1:consensusNumOfIter
            tmpA = selectNeighbor(G,0,ver);
            %P = 0.5*P*tmpA';
            s = 0.5*tmpA'*s;
            w = 0.5*tmpA'*w;
        end
        jointLikelihoodSums = s./w;
        %jointLikelihoodSums(1:end) = x;
        %jointLikelihoodSums = reshape(x,N,14);
        
        %%%%%%%%%%%%
        
end

end

function X = selectNeighbor(Y,type,matlabVersion)

Deg = sum(Y,2);
N = size(Y,1);
X = eye(N);

[a,dummy] = find(Y'==1);
switch type
    case 0
        ind = round(uniform(ones(1,N),Deg'));
    case 1
        %ind = round((Deg+1)'/2+Deg'+randn(1)*(Deg-1)');
        %ii = find( ind < 1);
        %ind(ii) = Deg(1);
        TODO
        
end

aa = cumsum([0,Deg']);
ind = ind+aa(1:end-1);
ind2 = (0:(N-1))*N+a(ind)';
 
X(ind2)=1;
X=X';

%
% for i=1:N
%     
%     a = find(Y(i,:)==1);
%     
%     switch type
%         case 0
%             if (matlabVersion == 1)
%                 [~,j] = max(rand(Deg(i),1));
%             else
%                 [dummy,j] = max(rand(Deg(i),1));
%                 
%             end
%         case 1
%             if (matlabVersion == 1)           
%                 [~,j] = max(randn(Deg(i),1));
%             else
%                 [dummy,j] = max(randn(Deg(i),1));
%             end
%     end
%     
%     X(i,a(j)) = 1;
%     
% end

end


