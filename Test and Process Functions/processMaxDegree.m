warning('off','all');
clear;clc;

filepath = 'Individual PF Results\LCpf_llh\';

% Number of particles for the filter
maxDegree_vector = [1:9];

% Number of random trials
no_trials = 200; 

RMSEFull = [];
steptimeFull = [];
weight_difFull = [];
N_effFull = [];
aggregate_error_ratioFull = [];

% xlabel parameters for the plot
xticklabel = {};
xtick = [];
groupSeparator = [];

colorgroup = [];

N = 500;

for i=1:numel(maxDegree_vector)
    % Load the tracking results
    filename = [filepath, 'Track3_LCpf'];
    filename = [filename,'_maxDegree',num2str(maxDegree_vector(i))];
    filename = [filename,'_N',num2str(N)];
    filename = [filename,'_trials',num2str(no_trials)];
    filename = [filename,'.mat']

    data = load(filename);
    [RMSEFull_sf, runtimeFull_sf, weight_dif_full_sf, N_effFull_sf, aggregate_error_ratio_sf, details{i}] = extractResults(data.results, data.parameters);

    RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
    steptimeFull = cat(4, steptimeFull, runtimeFull_sf);
    weight_difFull = cat(4, weight_difFull, mean(abs(weight_dif_full_sf),4));
    N_effFull = cat(4, N_effFull, N_effFull_sf);
    aggregate_error_ratioFull = cat(4, aggregate_error_ratioFull, aggregate_error_ratio_sf);
end

figure();
set(gcf,'color','white');
ylabel('RMSE');
xlabel('NGossip');
boxplot(squeeze((mean(weight_difFull,2))));
set(gca,'xtick', 1:numel(maxDegree_vector));
set(gca,'xticklabel', {'1', '2', '3', '4', '5', '6', '7', '8', '9'});
set(gca,'fontsize',45);
xlabel('max degree d');
ylabel('Relative error ratio')

% for i=1:numel(maxDegree_vector)
%     Hx_ss_dif(:,i) = mean(mean(mean(abs(details{i}.LCpf.Hx_ss_dif),4),3),1)';
% end
% figure();
% set(gcf,'color','white');
% boxplot(Hx_ss_dif/pi*180);%, 'colorgroup', colorgroup);
% ylabel({'Measurment approximation', ' error (degree)'});
% xlabel('d');
% set(gca,'xtick', 1:9);
% set(gca,'fontsize',32);

figure();
set(gcf,'color','white');
boxplot(squeeze(mean(RMSEFull,2)), 'colorgroup', colorgroup);
ylabel('RMSE');
xlabel('Max degree');
set(gca,'xtick', 1:9);
set(gca,'fontsize',32);
ylim([1,3.5]);

% figure();
% set(gcf,'color','white');
% boxplot(squeeze(sum(steptimeFull,2)), 'colorgroup', colorgroup);
% ylabel('Total runtime');
% xlabel('Max degree');
% set(gca,'xtick', 1:9);
% set(gca,'fontsize',32);