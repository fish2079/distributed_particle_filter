warning('off','all');
clear;clc;

path = '';

sigma_vector = [1, 2.5, 7.5, 10];

% Number of random trials
sim_parameters.no_trials = 80; 

RMSEFull = [];
runtimeFull = [];
weight_difFull = [];
N_effFull = [];
aggregate_error_ratioFull = [];

% xlabel parameters for the plot
xticklabel = {};
xtick = [];
groupSeparator = [];

colorgroup = [];

% Loop through each choice of particle number
for i=1:numel(sigma_vector)
    % Set number of particles
    sim_parameters.sigma = sigma_vector(i); 
    
    % Store the tracking results
    filename{i} = ['Track3_gossip_sigma',num2str(sim_parameters.sigma*10),'_trials',num2str(sim_parameters.no_trials),'.mat'];
    
    data = load(filename{i});
    [RMSEFull_sf, runtimeFull_sf, weight_dif_full_sf, N_effFull_sf, aggregate_error_ratio_sf, details{i}] = extractResults(data.results, data.parameters);
   
    RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
    runtimeFull = cat(4, runtimeFull, runtimeFull_sf);
    weight_difFull = cat(4, weight_difFull, mean(abs(weight_dif_full_sf),4));
    N_effFull = cat(4, N_effFull, N_effFull_sf);
    aggregate_error_ratioFull = cat(4, aggregate_error_ratioFull, aggregate_error_ratio_sf);
    
    xticklabel = [xticklabel, num2str(sigma_vector(i))];
    xtick = [xtick, size(RMSEFull_sf,1)*(i-1)+size(RMSEFull_sf,1)/2+0.5];
    groupSeparator(i) = size(RMSEFull_sf,1)*i+0.5;
    
    colorgroup = [colorgroup, 1:size(RMSEFull_sf,1)];
end

% figure();
% set(gcf,'color','white');
% boxplot(RMSE', 'colorgroup', colorgroup);
% ylabel('RMSE');
% xlabel('N');
% hold on;
% for i=1:numel(groupSeparator)-1
%     plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
% end
% set(gca,'xtick', xtick);
% set(gca,'xticklabel', xticklabel);
% set(gca,'fontsize',32);