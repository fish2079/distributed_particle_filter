function processClusterpf(track)
filepath = 'New Initialization Full Simulation Results\Clusterpf_param\'; 

if (track == 2)
    BSpf = 0.9304;
else
    BSpf = 1.6793;
end

k_vector = [4, 6, 9,16, 25, 36]; %[6, 9, 16, 25, 50];

% Number of gossip iterations
gossip_vector = [1,5, 10, 25, 50, 75, 100,200];

% Number of random trials
no_trials = 200; 

RMSEFull = [];

% xlabel parameters for the plot
xticklabel = {};

N = 1000;

% Loop through each choice of particle number
for i=1:numel(gossip_vector)
    for j=1:numel(k_vector)
        % Load the tracking results
        filename = [filepath, 'Track',num2str(track),'_Clusterpf'];
        filename = [filename, '_gossip',num2str(gossip_vector(i))];
        filename = [filename,'_K',num2str(k_vector(j))];
        filename = [filename,'_N',num2str(N)];
        filename = [filename,'_trials',num2str(no_trials)];
        filename = [filename,'.mat'];
        
        data = load(filename);
        [RMSEFull_sf] = extractResults(data.results, data.parameters, false);

        RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
    end
    xticklabel = [xticklabel, num2str(gossip_vector(i))];
end

[a,b]=meshgrid(gossip_vector, k_vector);
x_values = a.*b;

figure();
set(gcf,'color','white');
ylabel('RMSE');
xlabel('NGossip');
hold on;
RMSE_plot = reshape(squeeze(mean(mean(RMSEFull,2),3)), numel(k_vector), numel(gossip_vector));
legendText = {};
for i=1:numel(k_vector)
    plot(x_values(i,:), RMSE_plot(i,:), 'linewidth',6);
    legendText{i} = ['K = ', num2str(k_vector(i))];
end
plot(1:1000, BSpf*ones(1,1000),'k', 'linewidth',6);
legendText{i+1} = 'BSpf';
legend(legendText);
% set(gca,'xtick', 1:numel(gossip_vector));
% set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);
% xlim([1,numel(gossip_vector)]);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')
xlim([1,1000]);