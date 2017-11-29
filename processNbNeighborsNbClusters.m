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
KNN_vector = [3:10];
nbClusters_vector = 500;%[10, 20, 50, 100, 200, 500, 1000];

% Number of random trials
sim_parameters.no_trials = 200; 

RMSE = zeros(numel(KNN_vector), numel(nbClusters_vector));
time = zeros(numel(KNN_vector), numel(nbClusters_vector));

RMSEFull = [];
timeFull = [];
for j=1:numel(nbClusters_vector)
    sim_parameters.nbClusters = nbClusters_vector(j);
    xticklabel{j} = num2str(sim_parameters.nbClusters);
    for i=1:numel(KNN_vector)
        % Set number of particles
        sim_parameters.KNN = KNN_vector(i); 

        filename{1} = [filepath,'Clusterpf_KNN',num2str(sim_parameters.KNN),'_nbClusters',num2str(sim_parameters.nbClusters),'_trials',num2str(sim_parameters.no_trials),'.mat'];

        [RMSE_vector, time_vector] = processResults(filename, 1);
        RMSEFull = [RMSEFull; RMSE_vector];
        
        RMSE(i,j) = mean(RMSE_vector);
        time(i,j) = mean(time_vector);
        timeFull = [timeFull; time_vector];
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