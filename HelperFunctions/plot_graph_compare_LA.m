KNN_vector = [5,10,25,50,100,200,500,1000]; 
Epsilon_vector = [1/8, 1/6, 1/5, 1/4, 1/3, 1/2];

RMSE_plot = [RMSE_DT; RMSE_KNN(3:5,:)];% RMSE_Epsilon([1,4,6],:)];
figure();
set(gcf,'color','white');

% plot(1:12, 0.8238*ones(1,12),'k','linewidth',6);
plot(1:15, 1.5587*ones(1,15),'k','linewidth',6);
hold on;
plot(RMSE_plot', 'linewidth',6);
xlabel('m');
ylabel('RMSE');
legend('BSpf', 'DT', 'KNN-25','KNN-50','KNN-100'); %,'Proximity graph, \epsilon=1/8','Proximity graph, \epsilon=1/4','Proximity graph, \epsilon=1/2');
set(gca,'XTick',1:15);
set(gca,'XTicklabel',{'4','6','9','25','50','100','200','300','400','500','600','700','800','900','1000'});
set(gca,'fontsize',32);
% ylim([1.5,2]);
xlim([1,15]);

% axes('Position',[.7 .7 .2 .2])
% box on
% plot(1:3, 0.8238*ones(1,3),'k','linewidth',6);
% hold on;
% plot(RMSE_plot(1:4,1:3)','linewidth',6);
% set(gca,'fontsize',32);

weight_plot = [weight_DT; weight_KNN(3:5,:)];% weight_Epsilon([1,4,6],:)];
figure();
set(gcf,'color','white');
% plot(1:12, 1.5587*ones(1,12),'k','linewidth',6);
hold on;
plot(RMSE_plot', 'linewidth',6);
xlabel('m');
ylabel('Weight discrepancy');
legend('DT', 'KNN-25','KNN-50','KNN-100'); %,'Proximity graph, \epsilon=1/8','Proximity graph, \epsilon=1/4','Proximity graph, \epsilon=1/2');
set(gca,'XTick',1:15);
set(gca,'XTicklabel',{'4','6','9','25','50','100','200','300','400','500','600','700','800','900','1000'});
set(gca,'fontsize',32);
% ylim([1.5,2]);
xlim([1,15]);

% axes('Position',[.7 .7 .2 .2])
% box on
% plot(weight_plot(1:4,1:3)','linewidth',6);
% set(gca,'fontsize',32);