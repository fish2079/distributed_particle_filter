function processN(track)
path = 'New Initialization Full Simulation Results\BSpf_N\'; %'BSpf_Nb_particles_results\N_Track3_Range\';
N_vector = 100:100:2000;

% Number of random trials
sim_parameters.no_trials = 200; 

RMSEFull = [];
% steptimeFull = [];
% weight_difFull = [];
% N_effFull = [];
% aggregate_error_ratioFull = [];

% xlabel parameters for the plot
xticklabel = {};
xtick = [];
% groupSeparator = [];

% colorgroup = [];

% Loop through each choice of particle number
for i=1:numel(N_vector)
    % Set number of particles
    sim_parameters.N = N_vector(i); 
    
    % Store the tracking results
    filename{i} = [path, 'Track',num2str(track),'_BSpf_N',num2str(sim_parameters.N),'_trials',num2str(sim_parameters.no_trials),'.mat'];
    
    data = load(filename{i});
%     [RMSEFull_sf, runtimeFull_sf, weight_dif_full_sf, N_effFull_sf, aggregate_error_ratio_sf, details{i}] = extractResults(data.results, data.parameters);
    [RMSEFull_sf] = extractResults(data.results, data.parameters);
   
    RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
%     steptimeFull = cat(4, steptimeFull, runtimeFull_sf);
%     weight_difFull = cat(4, weight_difFull, mean(abs(weight_dif_full_sf),4));
%     N_effFull = cat(4, N_effFull, N_effFull_sf);
%     aggregate_error_ratioFull = cat(4, aggregate_error_ratioFull, aggregate_error_ratio_sf);
    
    xticklabel = [xticklabel, num2str(N_vector(i))];
    xtick = [xtick, size(RMSEFull_sf,1)*(i-1)+size(RMSEFull_sf,1)/2+0.5];
%     groupSeparator(i) = size(RMSEFull_sf,1)*i+0.5;
    
%     colorgroup = [colorgroup, 1:size(RMSEFull_sf,1)];
end

RMSE = squeeze(mean(mean(RMSEFull,3),2));
figure();
set(gcf,'color','white');
plot(RMSE,'c','linewidth',6);
% boxplot(RMSE', 'colorgroup', colorgroup);
ylabel('RMSE');
xlabel('N');
hold on;
% for i=1:numel(groupSeparator)-1
%     plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
% end
set(gca,'xtick', xtick(2:2:end));
set(gca,'xticklabel', xticklabel(2:2:end));
set(gca,'fontsize',32);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

RMSE_time = squeeze(mean(RMSEFull,3));
figure();
set(gcf,'color','white');
% plot(RMSE,'c','linewidth',6);
% boxplot(RMSE', 'colorgroup', colorgroup);
ylabel('RMSE');
xlabel('N');
% hold on;
plot(1:50, RMSE_time(:,[6,10,15,20]),'linewidth',6);
% set(gca,'xtick', xtick);
% set(gca,'xticklabel', xticklabel);
legend(xticklabel([6,10,15,20]));
xlabel('time');
ylabel('RMSE');
set(gca,'fontsize',32);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')
