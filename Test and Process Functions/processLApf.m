function processLApf(track)

filepath = 'New Initialization Full Simulation Results/LApf_param/'; %'LApf_parameter_results/';

if (track == 2)
    BSpf = 0.9304;
else
    BSpf = 1.6793;
end

m_vector = [4, 6, 9, 16, 25, 36];%,100,150,250];

% Number of gossip iterations
gossip_vector = [1, 5, 10, 25, 50, 75, 100, 200];

% Number of random trials
no_trials = 200; 

RMSEFull = [];

% xlabel parameters for the plot
xticklabel = {};

N = 1000;
% Loop through each choice of particle number
for i=1:numel(gossip_vector)
    for j=1:numel(m_vector)
        % Load the tracking results
        filename = [filepath, 'Track',num2str(track),'_LApf'];
        filename = [filename, '_gossip',num2str(gossip_vector(i))];
        filename = [filename,'_m',num2str(m_vector(j))];
        filename = [filename,'_N',num2str(N)];
        filename = [filename,'_trials',num2str(no_trials)];
        filename = [filename,'.mat'];
        
        data = load(filename);
        [RMSEFull_sf] = extractResults(data.results, data.parameters, false);

        RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
%         steptimeFull = cat(4, steptimeFull, runtimeFull_sf);
%         weight_difFull = cat(4, weight_difFull, mean(abs(weight_dif_full_sf),4));
%         N_effFull = cat(4, N_effFull, N_effFull_sf);
%         aggregate_error_ratioFull = cat(4, aggregate_error_ratioFull, aggregate_error_ratio_sf);
    end
    xticklabel = [xticklabel, num2str(gossip_vector(i))];
end

[a,b]=meshgrid(gossip_vector, m_vector);
x_values = a.*b;

figure();
set(gcf,'color','white');
ylabel('RMSE');
xlabel('NGossip');
hold on;
RMSE_plot = reshape(squeeze(mean(mean(RMSEFull,2),3)), numel(m_vector), numel(gossip_vector));
legendText = {};
for i=1:numel(m_vector)
    plot(x_values(i,:), RMSE_plot(i,:), 'linewidth',6);
    legendText{i} = ['m = ', num2str(m_vector(i))];
end
plot(1:1000, BSpf*ones(1,1000),'k', 'linewidth',6);
legendText{i+1} = 'BSpf';
legend(legendText);
% set(gca,'xtick', 1:numel(gossip_vector));
% set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k');
xlim([1,1000]);