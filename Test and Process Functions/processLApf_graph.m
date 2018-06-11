function processLApf_graph(track, plot_Epsilon)

filepath = 'New Initialization Full Simulation Results/LApf_graph/';

if (track == 2)
    BSpf = 0.9304;
else
    BSpf = 1.6793;
end

% m_vector = [4, 6, 9, 16, 25, 50];%, 100:100:700];
m_vector = [4, 6, 9, 16, 25, 50, 100, 400, 700, 1000];

%KNN_vector = [25,50,100]; %[10,25,50,100,200]; 
KNN_vector = [10, 25, 50, 100, 200]; 
Epsilon_vector = [1/8, 1/4, 1/2];

% Number of random trials
no_trials = 200; 

% Number of particles
N = 1000;

xticklabel = {};

% Extract Delaunay triangulation graph results
RMSEFull = [];
weight_difFull = [];
for j=1:numel(m_vector)
    % Load the tracking results
    filename = [filepath, 'Track',num2str(track),'_LApf_Weighted1'];
    filename = [filename, '_DT'];
    filename = [filename,'_m',num2str(m_vector(j))];
    filename = [filename,'_N',num2str(N)];
    filename = [filename,'_trials',num2str(no_trials)];
    filename = [filename,'.mat'];

    data = load(filename);
    [RMSEFull_sf, ~, weight_dif_full_sf] = extractResults(data.results, data.parameters, false);

    RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
    weight_difFull = cat(4, weight_difFull, mean(abs(weight_dif_full_sf),4));
    xticklabel = [xticklabel, num2str(m_vector(j))];
end

RMSEFull_DT = RMSEFull; 
weight_difFull_DT = weight_difFull;
RMSEFull = [];
weight_difFull = [];

% Extract KNN graph results
for j=1:numel(m_vector)
    for i=1:numel(KNN_vector)
        % Load the tracking results
        filename = [filepath, 'Track',num2str(track),'_LApf_Weighted1'];
        filename = [filename, '_KNN',num2str(KNN_vector(i))];
%         filename = [filename, '_Epsilon',num2str(round(1/Epsilon_vector(i)))];
%         filename = [filename, '_gossip',num2str(gossip_vector(i))];
        filename = [filename,'_m',num2str(m_vector(j))];
        filename = [filename,'_N',num2str(N)];
        filename = [filename,'_trials',num2str(no_trials)];
        filename = [filename,'.mat'];
        
        data = load(filename);
        [RMSEFull_sf, ~, weight_dif_full_sf] = extractResults(data.results, data.parameters, false);

        RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
%         steptimeFull = cat(4, steptimeFull, runtimeFull_sf);
        weight_difFull = cat(4, weight_difFull, mean(abs(weight_dif_full_sf),4));
%         N_effFull = cat(4, N_effFull, N_effFull_sf);
%         aggregate_error_ratioFull = cat(4, aggregate_error_ratioFull, aggregate_error_ratio_sf);
    end
end
RMSEFull_KNN = RMSEFull; 
weight_difFull_KNN = weight_difFull;
RMSEFull = [];
weight_difFull = [];

% Extract Epsilon graph results if necessary
if (plot_Epsilon)
    % Extract KNN graph results
    for j=1:numel(m_vector)
        for i=1:numel(Epsilon_vector)
            % Load the tracking results
            filename = [filepath, 'Track',num2str(track),'_LApf_Weighted1'];
            filename = [filename, '_Epsilon',num2str(round(1/Epsilon_vector(i)))];
            filename = [filename,'_m',num2str(m_vector(j))];
            filename = [filename,'_N',num2str(N)];
            filename = [filename,'_trials',num2str(no_trials)];
            filename = [filename,'.mat'];

            data = load(filename);
            [RMSEFull_sf, ~, weight_dif_full_sf] = extractResults(data.results, data.parameters, false);

            RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
            weight_difFull = cat(4, weight_difFull, mean(abs(weight_dif_full_sf),4));
        end
    end
    RMSEFull_Epsilon = RMSEFull; 
    weight_difFull_Epsilon = weight_difFull;
    RMSEFull = [];
    weight_difFull = [];
end

% Plot the figure of RMSE
figure();
set(gcf,'color','white');
ylabel('RMSE');
xlabel('m');
hold on;
RMSE_plot_DT = squeeze(mean(mean(RMSEFull_DT,2),3))';
RMSE_plot_KNN = reshape(squeeze(mean(mean(RMSEFull_KNN,2),3)), numel(KNN_vector), numel(m_vector));
RMSE_plot = [RMSE_plot_DT; RMSE_plot_KNN];
if (plot_Epsilon)
    RMSE_plot_Epsilon = reshape(squeeze(mean(mean(RMSEFull_Epsilon,2),3)), numel(Epsilon_vector), numel(m_vector));
    RMSE_plot = [RMSE_plot; RMSE_plot_Epsilon];
end
plot(m_vector, RMSE_plot_DT,'k','linewidth',6);
legendText{1} = 'Delaunay';
hold on;
for i=1:numel(KNN_vector)
    plot(m_vector, RMSE_plot_KNN(i,:),'linewidth',6);
    legendText{i+1} = ['KNN-', num2str(KNN_vector(i))];
end

if (plot_Epsilon)
    for i=1:numel(Epsilon_vector)
        plot(m_vector, RMSE_plot_Epsilon(i,:),'linewidth',6);
        legendText{i+1+numel(KNN_vector)} = ['Proximity graph \epsilon=1/', num2str(1/Epsilon_vector(i))];
    end
end

legend(legendText);
% set(gca,'xtick', 1:numel(m_vector));
% set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);
xlim([4,50]);
% xlim([1,numel(m_vector)]);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

% Plot the figure of weight discrepancy
figure();
set(gcf,'color','white');
ylabel('Weight Discrepancy');
xlabel('m');
hold on;
weight_plot_DT = squeeze(mean(mean(weight_difFull_DT,2),3))';
weight_plot_KNN = reshape(squeeze(mean(mean(weight_difFull_KNN,2),3)), numel(KNN_vector), numel(m_vector));
if (plot_Epsilon)
    weight_plot_Epsilon = reshape(squeeze(mean(mean(weight_difFull_Epsilon,2),3)), numel(Epsilon_vector), numel(m_vector));
end
plot(m_vector, weight_plot_DT,'k','linewidth',6);
legendText{1} = 'Delaunay';
hold on;
for i=1:numel(KNN_vector)
    plot(m_vector, weight_plot_KNN(i,:),'linewidth',6);
    legendText{i+1} = ['KNN-', num2str(KNN_vector(i))];
end

if (plot_Epsilon)
    for i=1:numel(Epsilon_vector)
        plot(m_vector, weight_plot_Epsilon(i,:),'linewidth',6);
        legendText{i+1+numel(KNN_vector)} = ['Proximity graph \epsilon=1/', num2str(1/Epsilon_vector(i))];
    end
end

legend(legendText);
% set(gca,'xtick', 1:numel(m_vector));
% set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);
% xlim([1,numel(m_vector)]);
xlim([4,50]);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')