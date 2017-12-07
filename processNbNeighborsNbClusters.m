% Test script for comparing the algorithms' performance in terms of number 
% of particles on the simulated track
% For each parameter value, the script generates the simulated track,
% measurements and runs the specified tracking algorithms for a number of
% MC trials
% Note that all algorithms can only track one single target
% Multiple targets with target birth/death, clutter measurements are not
% considered

% Jun Ye Yu
% McGill University
% jun.y.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

warning('off','all');
clear;clc;

filepath = 'Results_Clusterpf\';
% filepath = '';

% Number of particles for the filter
KNN_vector = [3:10,20,50,100];
nbClusters_vector = [6, 10, 20, 50, 100, 200, 500, 1000];

% Number of random trials
sim_parameters.no_trials = 100; 

RMSE = zeros(numel(KNN_vector), numel(nbClusters_vector));
time = zeros(numel(KNN_vector), numel(nbClusters_vector));

RMSEFull = [];
timeFull = [];
KNNtimeFull = [];
clustertimeFull = [];
gammatimeFull = [];
steptimeFull = [];

for j=1:numel(nbClusters_vector)
    sim_parameters.nbClusters = nbClusters_vector(j);
    xticklabel{j} = num2str(sim_parameters.nbClusters);
    for i=1:numel(KNN_vector)
        % Set number of particles
        sim_parameters.KNN = KNN_vector(i); 

        filename{i} = [filepath,'Clusterpf_KNN',num2str(sim_parameters.KNN),'_nbClusters',num2str(sim_parameters.nbClusters),'_trials',num2str(sim_parameters.no_trials),'.mat'];

        [RMSE_vector, ~, time_vector, runtimeFull, detail] = extractResults(filename{i});
        
        RMSE(i,j) = nanmean(RMSE_vector);
        RMSEFull = [RMSEFull; RMSE_vector];
        
        timeFull = [timeFull; squeeze(sum(runtimeFull,2))'];
        time(i,j) = mean(squeeze(sum(runtimeFull,2))');
        
        KNNtimeFull = [KNNtimeFull; sum(detail.Clusterpf.KNN_time,1)];
        clustertimeFull = [clustertimeFull; sum(detail.Clusterpf.cluster_time,1)];
        gammatimeFull = [gammatimeFull; sum(detail.Clusterpf.gamma_time,1)];
        steptimeFull = [steptimeFull;squeeze(sum(runtimeFull,2))'];
    end
end

figure();
boxplot(RMSEFull');

figure();
boxplot(timeFull');
% figure();
% set(gcf,'color','white');
% boxplot(RMSEFull(1:4:end,:)');
% xlabel('Number of Clusters');
% ylabel('RMSE');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% ylim([0,8]);
% 
% figure();
% set(gcf,'color','white');
% boxplot(RMSEFull(2:4:end,:)');
% xlabel('Number of Clusters');
% ylabel('RMSE');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% ylim([0,8]);
% 
% figure();
% set(gcf,'color','white');
% boxplot(RMSEFull(3:4:end,:)');
% xlabel('Number of Clusters');
% ylabel('RMSE');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% ylim([0,8]);
% 
% figure();
% set(gcf,'color','white');
% boxplot(timeFull(1:4:end,:)');
% xlabel('Number of Clusters');
% ylabel('Total runtime');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% 
% figure();
% set(gcf,'color','white');
% boxplot(timeFull(2:4:end,:)');
% xlabel('Number of Clusters');
% ylabel('Total runtime');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% 
% figure();
% set(gcf,'color','white');
% boxplot(timeFull(3:4:end,:)');
% xlabel('Number of Clusters');
% ylabel('Total runtime');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));