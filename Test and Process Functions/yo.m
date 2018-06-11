figure();
set(gcf,'color','white');
xlabel('Scalars transmitted per sensor');
ylabel('RMSE')
hold on;

gossip_vector = [1,5,10,25,50,75,100,200];

% Plot CSSpf
plot(x_values_CSSpf, RMSE_CSS,'k-+','linewidth',6,'markersize',20);
legendText = {'CSS[5] (6)'};
legendIndex = 2;

% Plot GApf
plot(x_values_GApf, RMSE_GA, 'g-+','linewidth',6,'markersize',20);
legendText{legendIndex} = ['GApf[3] (14)'];
legendIndex = legendIndex+1;

% Plot LCpf
% set(gca,'ColorOrderIndex', 1);
scalars = [4,9,16,25,36];
for i=1%:size(x_values_LCpf,1)
    plot(x_values_LCpf(i,:), RMSE_LCpf_track3(i,:),'-+','linewidth',6,'markersize',20);
    legendText{legendIndex} = ['LCpf[4] (',num2str(scalars(i)),')'];
    legendIndex = legendIndex+1;
end

% Plot LCpf-GS
% set(gca,'ColorOrderIndex', 1);
scalars = [4,9,16,25,36];
for i=1%1:size(x_values_LCpf,1)
    plot(x_values_LCpf(i,:), RMSE_LCpf_GS_track3(i,:),'-+','linewidth',6,'markersize',20);
    legendText{legendIndex} = ['LCpf-GS (',num2str(scalars(i)),')'];
    legendIndex = legendIndex+1;
end

% Plot LApf
% set(gca,'ColorOrderIndex', 1);
scalars = [4,6,9,16,25,36];
for i=3%1:size(x_values_LApf,1)
    plot(x_values_LApf(i,:), RMSE_LApf_track3(i,:),'-+','linewidth',6,'markersize',20);
    legendText{legendIndex} = ['LApf (',num2str(scalars(i)),')'];
    legendIndex = legendIndex+1;
end

% Plot Clusterpf
% set(gca,'ColorOrderIndex', 1);
scalars = [4,6,9,16,25,36];
for i=3%1:size(x_values_Clusterpf,1)
    plot(x_values_Clusterpf(i,:), RMSE_Clusterpf_track3(i,:),'-+','linewidth',6,'markersize',20);
    legendText{legendIndex} = ['Clusterpf (',num2str(scalars(i)),')'];
    legendIndex = legendIndex+1;
end

legend(legendText);
xlim([1,500]);
set(gca,'fontsize',32);

BSpf = 1.6793; %0.9304;
plot(1:500, ones(1,500)*BSpf,'m','linewidth',6);

grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

set(gca,'YScale','log');