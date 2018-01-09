warning('off','all');
clear;clc;

filepath = 'Individual PF Results\LCpf\';

% Number of particles for the filter
N_vector = [100, 250, 500, 1000];

% Number of gossip iterations
gossip_vector = [10, 25, 50, 100, 200];


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
max_degree = 1;

% Loop through each choice of particle number
for i=1:numel(gossip_vector)
    for j=1:numel(N_vector)
        % Load the tracking results
        filename = [filepath, 'Track3_CSSpf'];
        filename = [filename, '_gossip',num2str(gossip_vector(i))];
        filename = [filename,'_maxDegree',num2str(max_degree)];
        filename = [filename,'_N',num2str(N_vector(j))];
        filename = [filename,'_trials',num2str(no_trials)];
        filename = [filename,'.mat'];
        
        data = load(filename);
        [RMSEFull_sf, runtimeFull_sf, weight_dif_full_sf, N_effFull_sf, aggregate_error_ratio_sf, details{i}] = extractResults(data.results, data.parameters);

        RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
        steptimeFull = cat(4, steptimeFull, runtimeFull_sf);
        weight_difFull = cat(4, weight_difFull, mean(abs(weight_dif_full_sf),4));
        N_effFull = cat(4, N_effFull, N_effFull_sf);
        aggregate_error_ratioFull = cat(4, aggregate_error_ratioFull, aggregate_error_ratio_sf);
    end
    xticklabel = [xticklabel, num2str(gossip_vector(i))];
    xtick = [xtick, numel(N_vector)*(i-1)+numel(N_vector)/2+0.5];
    groupSeparator(i) = numel(N_vector)*i+0.5;
    
    colorgroup = [colorgroup, 1:numel(N_vector)];
end

figure();
set(gcf,'color','white');
boxplot(squeeze(mean(RMSEFull,2)), 'colorgroup', colorgroup);
ylabel('RMSE');
xlabel('NGossip');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,90],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);
ylim([0,90]);

figure();
set(gcf,'color','white');
boxplot(squeeze(sum(steptimeFull,2)), 'colorgroup', colorgroup);
ylabel('Total runtime');
xlabel('NGossip');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);

figure();
set(gcf,'color','white');
boxplot(squeeze(sum(aggregate_error_ratioFull,2)), 'colorgroup', colorgroup);
ylabel('Aggregate error ratio');
xlabel('NGossip');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,10^4],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);
