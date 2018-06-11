figure();
set(gcf,'color','white');

plot(1:8, 0.8238*ones(1,8),'k','linewidth',6);
% plot(1:12, 1.5587*ones(1,12),'k','linewidth',6);
hold on;
plot(RMSE(:,2:5),'linewidth',6);
% plot(RMSE(:,1),'m','linewidth',6);
% plot(RMSE(:,2),'r','linewidth',6);
% plot(RMSE(:,3),'g','linewidth',6);
% plot(RMSE(:,4),'b','linewidth',6);
% plot(RMSE(:,5),'m-o','linewidth',6,'markersize',20);
% plot(RMSE(:,6),'r-o','linewidth',6,'markersize',20);
% plot(RMSE(:,7),'g-o','linewidth',6,'markersize',20);
xlabel('NGossip');
ylabel('RMSE');
% legend('BSpf', 'K=6', 'K=9', 'K=25', 'K=50', 'K=100', 'K=250', 'K=300');
legend('BSpf', 'K=9', 'K=25', 'K=50', 'K=100');
set(gca,'XTick',1:8);
set(gca,'XTicklabel',{'1','5', '10', '25', '50', '75', '100', '200'});
set(gca,'fontsize',32);
% ylim([1.5,2]);
xlim([1,8]);