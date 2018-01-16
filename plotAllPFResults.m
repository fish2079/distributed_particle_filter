clear;clc;
load('All PF Results\Track2\Track2_results2.mat');

linewidth = 6;
nb_data = size(RMSE,2);

figure();
set(gcf,'color','white');
plot(1:nb_data, RMSE(1,:), 'm', 'linewidth',linewidth);
hold on;
plot(1:nb_data, RMSE(2,:), 'g', 'linewidth',linewidth);
plot(1:nb_data, RMSE(3,:), 'rd', 'linewidth',linewidth);
plot(1:nb_data, RMSE(4,:), 'bd', 'linewidth',linewidth);
plot(1:nb_data, RMSE(5,:), 'r', 'linewidth',linewidth);
plot(1:nb_data, RMSE(6,:), 'b', 'linewidth',linewidth);
plot(1:nb_data, ones(1,nb_data)*BSpf_RMSE, 'k', 'linewidth',linewidth);
% ylim([1.5,3]);
% xlim([1,20]);
set(gca, 'xtick',1:2:20);
set(gca, 'xticklabel', [60:120:1200]);
legend('CSSpf', 'LCpf', 'LApf (m=6)', 'Clusterpf (C=6)', 'LApf (NGossip=60)', 'Clusterpf (NGossip=60)','BSpf');
xlabel('Total communication overhead per time step');
ylabel('RMSE');
set(gca,'fontsize',32);

figure();
set(gcf,'color','white');
plot(1:nb_data, time(1,:), 'm', 'linewidth',linewidth);
hold on;
plot(1:nb_data, time(2,:), 'g', 'linewidth',linewidth);
plot(1:nb_data, time(3,:), 'rd', 'linewidth',linewidth);
plot(1:nb_data, time(4,:), 'bd', 'linewidth',linewidth);
plot(1:nb_data, time(5,:), 'r', 'linewidth',linewidth);
plot(1:nb_data, time(6,:), 'b', 'linewidth',linewidth);
plot(1:nb_data, ones(1,nb_data)*BSpf_time, 'k', 'linewidth',linewidth);
% ylim([1.5,3]);
% xlim([1,20]);
set(gca, 'xtick',1:2:20);
set(gca, 'xticklabel', [60:120:1200]);
legend('CSSpf', 'LCpf', 'LApf (m=6)', 'Clusterpf (C=6)', 'LApf (NGossip=60)', 'Clusterpf (NGossip=60)','BSpf');
xlabel('Total communication overhead per time step');
ylabel('Total runtime (s)');
set(gca,'fontsize',32);