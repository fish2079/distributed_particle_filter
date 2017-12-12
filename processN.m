warning('off','all');
clear;clc;

path = 'Results_N\Track3\';
% path = '';
% Number of particles for the filter
N_vector = [500, 750, 1000, 1250, 1500];

% Number of random trials
sim_parameters.no_trials = 500; 

RMSE = [];
RMSEFull = [];
runtime = [];
runtimeFull = [];

% xlabel parameters for the plot
xticklabel = {};
xtick = [];
groupSeparator = [];

colorgroup = [];

% Loop through each choice of particle number
for i=1:numel(N_vector)
    % Set number of particles
    sim_parameters.N = N_vector(i); 
    
    % Store the tracking results
    filename{i} = [path, 'Track3_N',num2str(sim_parameters.N),'_trials',num2str(sim_parameters.no_trials),'.mat'];
    
    [RMSE_sf, RMSEFull_sf, runtime_sf, runtimeFull_sf, detail{i}] = extractResults(filename{i});
    
    RMSE = [RMSE; RMSE_sf];
    RMSEFull = [RMSEFull; RMSEFull_sf];
    runtime = [runtime; runtime_sf];
    runtimeFull = [runtimeFull; runtimeFull_sf];
    
    xticklabel = [xticklabel, num2str(N_vector(i))];
    xtick = [xtick, size(RMSE_sf,1)*(i-1)+size(RMSE_sf,1)/2+0.5];
    groupSeparator(i) = size(RMSE_sf,1)*i+0.5;
    
    colorgroup = [colorgroup, 1:size(RMSE_sf,1)];
end

figure();
set(gcf,'color','white');
boxplot(RMSE', 'colorgroup', colorgroup);
ylabel('RMSE');
xlabel('N');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);

% figure();
% boxplot(runtime', 'colorgroup', colorgroup);
% ylabel('total runtime');
% xlabel('N');
% hold on;
% for i=1:numel(groupSeparator)-1
%     plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
% end
% set(gca,'xtick', xtick);
% set(gca,'xticklabel', xticklabel);
% set(gca,'fontsize',32);