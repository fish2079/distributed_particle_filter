KNN_vector = [5,10,25,50,100,200,500,1000]; 

figure();
set(gcf,'color','white');
RMSE_plot = [RMSE_DT; RMSE_KNN(4:6,:)];% RMSE_Epsilon([1,4,6],:)];

% plot(1:12, 0.8238*ones(1,12),'k','linewidth',6);
plot(1:15, 1.5587*ones(1,15),'k','linewidth',6);
hold on;
plot(RMSE_plot', 'linewidth',6);
xlabel('K');
ylabel('RMSE');
legend('BSpf', 'DT', 'KNN-50','KNN-100', 'KNN-200'); %,'Proximity graph, \epsilon=1/8','Proximity graph, \epsilon=1/4','Proximity graph, \epsilon=1/2');
set(gca,'XTick',1:12);
set(gca,'XTicklabel',{'4','6','9','15','25','30','50','100','150','200','250','300'});
set(gca,'fontsize',32);
% ylim([1.5,2]);
xlim([1,12]);

figure();
set(gcf,'color','white');
weight_plot = [weight_DT; weight_KNN(4:6,:)];% RMSE_Epsilon([1,4,6],:)];

hold on;
plot(weight_plot', 'linewidth',6);
xlabel('K');
ylabel('Weight discprenacy');
legend('DT', 'KNN-50','KNN-100', 'KNN-200'); %,'Proximity graph, \epsilon=1/8','Proximity graph, \epsilon=1/4','Proximity graph, \epsilon=1/2');
set(gca,'XTick',1:12);
set(gca,'XTicklabel',{'4','6','9','15','25','30','50','100','150','200','250','300'});
set(gca,'fontsize',32);
% ylim([1.5,2]);
xlim([1,12]);