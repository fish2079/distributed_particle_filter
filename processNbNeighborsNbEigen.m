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

filepath = 'OldResults\Results_LApf\';
% filepath = '';
% Number of particles for the filter
KNN_vector = [10];
nbEig_vector = [5, 10, 15, 20, 50, 100, 250, 500, 750];

% Number of random trials
sim_parameters.no_trials = 200; 

RMSE = zeros(numel(KNN_vector), numel(nbEig_vector));
time = zeros(numel(KNN_vector), numel(nbEig_vector));
RMSEFull = [];

RMSETime = [];
timeFull = [];
for j=1:numel(nbEig_vector)
    sim_parameters.nbEig = nbEig_vector(j);
    xticklabel{j} = num2str(sim_parameters.nbEig);
    for i=1:numel(KNN_vector)
        % Set number of particles
        sim_parameters.KNN = KNN_vector(i); 

        filename{1} = [filepath,'KNN',num2str(sim_parameters.KNN),'_nbEig',num2str(sim_parameters.nbEig),'_trials',num2str(sim_parameters.no_trials),'.mat'];

        [RMSE_vector, time_vector] = processResults(filename, 1);
        
        RMSE(i,j) = nanmean(RMSE_vector);
        RMSEFull = [RMSEFull; RMSE_vector];
        RMSESTD(i,j) = std(RMSE_vector);
        
        timeFull = [timeFull; time_vector];
        time(i,j) = mean(time_vector);
        timeSTD(i,j) = std(time_vector);
    end
end

% figure();
% set(gcf,'color','white');
% boxplot(RMSEFull(1:4:end,:)');
% xlabel('Number of Eigenvectors');
% ylabel('RMSE');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% ylim([0,8]);
% 
% figure();
% set(gcf,'color','white');
% boxplot(RMSEFull(2:4:end,:)');
% xlabel('Number of Eigenvectors');
% ylabel('RMSE');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% ylim([0,8]);
% 
% figure();
% set(gcf,'color','white');
% boxplot(RMSEFull(3:4:end,:)');
% xlabel('Number of Eigenvectors');
% ylabel('RMSE');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% ylim([0,8]);
% 
% figure();
% set(gcf,'color','white');
% boxplot(timeFull(1:4:end,:)');
% xlabel('Number of Eigenvectors');
% ylabel('Total runtime');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% 
% figure();
% set(gcf,'color','white');
% boxplot(timeFull(2:4:end,:)');
% xlabel('Number of Eigenvectors');
% ylabel('Total runtime');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));
% 
% figure();
% set(gcf,'color','white');
% boxplot(timeFull(3:4:end,:)');
% xlabel('Number of Eigenvectors');
% ylabel('Total runtime');
% set(gca,'fontsize',48);
% set(gca,'xtick',2:2:numel(xticklabel));
% set(gca,'xticklabel', xticklabel(2:2:end));