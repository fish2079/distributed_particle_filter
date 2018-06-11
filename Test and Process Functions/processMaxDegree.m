function processMaxDegree(track, GS_flag)
filepath = 'New Initialization Full Simulation Results/LCpf_param/'; %'Individual PF Results\LCpf_llh\';

% Number of particles for the filter
maxDegree_vector = [1,2,3,4,5];

gossip_vector = [1,5,10,25,50,75,100,200];

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

N = 1000;

for j=1:numel(gossip_vector)
    for i=1:numel(maxDegree_vector)
        % Load the tracking results
        if (GS_flag)
            filename = [filepath, 'Track',num2str(track),'_LCpf-GS_'];
        else
            filename = [filepath, 'Track',num2str(track),'_LCpf'];
        end
        filename = [filename,'_gossip',num2str(gossip_vector(j))];
        filename = [filename,'_maxDegree',num2str(maxDegree_vector(i))];
        filename = [filename,'_N',num2str(N)];
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
end

RMSE_plot = reshape(squeeze(mean(mean(RMSEFull,2),3)), numel(maxDegree_vector), numel(gossip_vector));

[a,b]=meshgrid(gossip_vector, maxDegree_vector);
x_values = a.*b;

figure();
hold on;
set(gcf,'color','white');
ylabel('RMSE');
xlabel('Scalars transmitted per sensor');
for i=1:numel(maxDegree_vector)
    plot(x_values(i,:), RMSE_plot(i,:), '+-','linewidth',6,'markersize',24);
end
legend('Nb of coefficients: 4','Nb of coefficients: 9','Nb of coefficients: 16','Nb of coefficients: 25','Nb of coefficients: 36');
set(gca,'fontsize',32);