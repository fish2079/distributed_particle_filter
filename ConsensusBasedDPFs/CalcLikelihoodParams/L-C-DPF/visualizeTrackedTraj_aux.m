load sensorsXY;
load adjMtx;
load x_n_viz;
load x_n_hat_viz;

numTargets = size(x_n_hat_viz,2)/4; %CURRENTLY THIS SCRIPT WORKS ONLY FOR UP TO 4 TARGETS !!!!!!!!!!!
simTime = size(x_n_hat_viz,3);
visualizedSensor = 1;

numSensors = size(sensorsXY,1);

x_n_hat = reshape(x_n_hat_viz(visualizedSensor,:,:),4*numTargets,simTime);

figure(59);
hold on;
plot(sensorsXY(:,1),sensorsXY(:,2),'o','LineWidth',8); 
for iter_sensor=1:numSensors 
    
    for iter_sensor2=1:numSensors 
        
        if adjMtx(iter_sensor,iter_sensor2) == 1
            
            line ([sensorsXY(iter_sensor,1) sensorsXY(iter_sensor2,1)],[sensorsXY(iter_sensor,2) sensorsXY(iter_sensor2,2)],'LineWidth',2,'LineStyle',':') ; 
            
        end
    end
end

%Target 1:
plot(x_n_viz(1,:),x_n_viz(2,:),'-','Color',[0.05,0.7,0.05],'LineWidth',4);
plot(x_n_hat(1,:),x_n_hat(2,:),'r-','LineWidth',4);

if numTargets > 1
%Target 2:
plot(x_n_viz(5,:),x_n_viz(6,:),'c-','LineWidth',4);
plot(x_n_hat(5,:),x_n_hat(6,:),'y-','LineWidth',4);
end

if numTargets > 2
%Target 3:
plot(x_n_viz(9,:),x_n_viz(10,:),'b-','LineWidth',4);
plot(x_n_hat(9,:),x_n_hat(10,:),'m-','LineWidth',4);
end

if numTargets > 3
%Target 4:
plot(x_n_viz(13,:),x_n_viz(14,:),'k-','LineWidth',4);
plot(x_n_hat(13,:),x_n_hat(14,:),'g-','LineWidth',4);
end

set(gca,'XColor','w','YColor','w');