clear;clc;close all;

load('testLApfGraphData.mat')

F.LA.KNN = 25;
F.LA.epsilon = 1/4;

d = size(x_predicted,1)-1;

% First have each sensor compute local log-likelihood using only local
% measurements
for i=1:numel(D.sensorID)   
    z_received = D.measurements(:,i);
    % Compute expected measurement
    z_expected = obs.model(x_predicted(1:d,:), D.sensorLoc(:,i), obs);
    
    % Compute the Gaussian log-likelihood
    z_dif = F.minus(z_received, z_expected);
    
    log_lh_ss(i,:) = log(mvnpdf(z_dif', obs.mu', obs.R)+realmin)';
end

% KNN
idx = knnsearch(x_predicted(1:2,:)', x_predicted(1:2,:)','k',F.LA.KNN+1);
idx = idx(:,2:end);
% Now construct the adjacency matrix
A_KNN = zeros(F.N, F.N);
for i=1:F.N
    % particle i is connected to its K nearest neighbor
    A_KNN(i,idx(i,:)) = 1;
    % the connection is symmetric
    A_KNN(idx(i,:),i) = 1;
end
        
A_Delaunay = DelaunayGraph(x_predicted(1:2,:)');

A_Epsilon = EpsilonGraph(x_predicted(1:2,:)', F.LA.epsilon);

A_Dummy = zeros(F.N, F.N);
for i=1:F.N-1
    A_Dummy(i,i+1) = 1;
    A_Dummy(i+1,i) = 1;
end
A_Dummy(F.N,1) = 1;
A_Dummy(1,F.N) = 1;


% figure();
% plot(graph(A_Delaunay),'XData',x_predicted(1,:), 'YData',x_predicted(2,:));
% 
% figure();
% plot(graph(A_KNN),'XData',x_predicted(1,:), 'YData',x_predicted(2,:));
% 
% figure();
% plot(graph(A_Epsilon),'XData',x_predicted(1,:), 'YData',x_predicted(2,:));
% 
% figure();
% plot(graph(A_Dummy),'XData',x_predicted(1,:), 'YData',x_predicted(2,:));

for i=1:F.N
    x1 = x_predicted(1:2,i);
    x2 = x_predicted(1:2,A_Delaunay(i,:)>0);
    dist = sqrt((x1(1,:)-x2(1,:)).^2+(x1(2,:)-x2(2,:)).^2);
    A_Delaunay(i,A_Delaunay(i,:)>0) = 1./dist;
    
    x1 = x_predicted(1:2,i);
    x2 = x_predicted(1:2,A_KNN(i,:)>0);
    dist = sqrt((x1(1,:)-x2(1,:)).^2+(x1(2,:)-x2(2,:)).^2);
    A_KNN(i,A_KNN(i,:)>0) = 1./dist;
    
    x1 = x_predicted(1:2,i);
    x2 = x_predicted(1:2,A_Epsilon(i,:)>0);
    dist = sqrt((x1(1,:)-x2(1,:)).^2+(x1(2,:)-x2(2,:)).^2);
    A_Epsilon(i,A_Epsilon(i,:)>0) = 1./dist;
    
    x1 = x_predicted(1:2,i);
    x2 = x_predicted(1:2,A_Dummy(i,:)>0);
    dist = sqrt((x1(1,:)-x2(1,:)).^2+(x1(2,:)-x2(2,:)).^2);
    A_Dummy(i,A_Dummy(i,:)>0) = 1./dist;
end

L_Delaunay = diag(sum(A_Delaunay,2)) - A_Delaunay;
L_KNN = diag(sum(A_KNN,2)) - A_KNN;
L_Epsilon = diag(sum(A_Epsilon,2)) - A_Epsilon;
L_Dummy = diag(sum(A_Dummy,2)) - A_Dummy;

[V_full_Delaunay,eig_Delaunay] = eig(L_Delaunay);
[V_full_KNN,eig_KNN] = eig(L_KNN);
[V_full_Epsilon,eig_Epsilon] = eig(L_Epsilon);
[V_full_Dummy,eig_Dummy] = eig(L_Dummy);

gamma_true = sum(log_lh_ss,1)';
gamma_true = gamma_true-max(gamma_true);
weight_true = exp(gamma_true);
weight_true = weight_true/sum(weight_true);

idx = 1;
for m=10:10:F.N
    V_Delaunay = V_full_Delaunay(:,1:m);
    V_KNN = V_full_KNN(:,1:m);
    V_Epsilon = V_full_Epsilon(:,1:m);
    V_Dummy = V_full_Dummy(:,1:m);
    
    alpha_Delaunay = sum(V_Delaunay'*log_lh_ss',2);
    alpha_KNN = sum(V_KNN'*log_lh_ss',2);
    alpha_Epsilon = sum(V_Epsilon'*log_lh_ss',2);
    alpha_Dummy = sum(V_Dummy'*log_lh_ss',2);
    
    std_alpha(1,idx) = norm((eye(m)-ones(m,m)/m)*sum(alpha_Delaunay,2))^2;
    std_alpha(2,idx) = norm((eye(m)-ones(m,m)/m)*sum(alpha_KNN,2))^2;
    std_alpha(3,idx) = norm((eye(m)-ones(m,m)/m)*sum(alpha_Epsilon,2))^2;
    std_alpha(4,idx) = norm((eye(m)-ones(m,m)/m)*sum(alpha_Dummy,2))^2;
        
    gamma_Delaunay = V_Delaunay*alpha_Delaunay;
    gamma_KNN = V_KNN*alpha_KNN;
    gamma_Epsilon = V_Epsilon*alpha_Epsilon;
    gamma_Dummy = V_Dummy*alpha_Dummy;
    
    gamma_Delaunay = gamma_Delaunay-max(gamma_Delaunay);
    weight_Delaunay = exp(gamma_Delaunay);
    weight_Delaunay = weight_Delaunay/sum(weight_Delaunay);
    
    gamma_KNN = gamma_KNN-max(gamma_KNN);
    weight_KNN = exp(gamma_KNN);
    weight_KNN = weight_KNN/sum(weight_KNN);
    
    gamma_Dummy = gamma_Dummy-max(gamma_Dummy);
    weight_Dummy = exp(gamma_Dummy);
    weight_Dummy = weight_Dummy/sum(weight_Dummy);
    
    gamma_Epsilon = gamma_Epsilon-max(gamma_Epsilon);
    weight_Epsilon = exp(gamma_Epsilon);
    weight_Epsilon = weight_Epsilon/sum(weight_Epsilon);
    weight_Dummy = weight_Dummy/sum(weight_Dummy);
    
    weightDif(:,idx) = [sum(abs(weight_true-weight_Delaunay)),sum(abs(weight_true-weight_KNN)),sum(abs(weight_true-weight_Epsilon)),sum(abs(weight_true-weight_Dummy))]';
    idx = idx+1;
end

figure();
plot(weightDif','linewidth',6);
legend('DT','KNN','Epsilon','Dummy');



