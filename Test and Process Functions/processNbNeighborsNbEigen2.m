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

filepath = 'Results_LApf\';

% Number of particles for the filter
KNN_vector = 3:10;
nbEig_vector = [6,10, 20, 50, 100, 200, 500, 1000];

% Number of random trials
sim_parameters.no_trials = 100; 

RMSE = zeros(numel(KNN_vector), numel(nbEig_vector));
time = zeros(numel(KNN_vector), numel(nbEig_vector));

RMSEFull = [];
timeFull = [];

KNNtimeFull = [];
eigtimeFull = [];
steptimeFull = [];

% xlabel parameters for the plot
xticklabel = {};
xtick = [];
groupSeparator = [];

for i=1:numel(KNN_vector)
    % Set number of neighbors
    sim_parameters.KNN = KNN_vector(i); 
        
    for j=1:numel(nbEig_vector)
        sim_parameters.nbEig = nbEig_vector(j);
        
        filename{i} = [filepath,'LApf_KNN',num2str(sim_parameters.KNN),'_nbEig',num2str(sim_parameters.nbEig),'_trials',num2str(sim_parameters.no_trials),'.mat'];

        [RMSE_vector, ~, time_vector, runtimeFull, detail] = extractResults(filename{i});
        
        RMSE(i,j) = nanmean(RMSE_vector);
        RMSEFull = [RMSEFull; RMSE_vector];
        
        
        timeFull = [timeFull; squeeze(sum(runtimeFull,2))'];
        time(i,j) = mean(squeeze(sum(runtimeFull,2))');
        
        KNNtimeFull = [KNNtimeFull; mean(detail.LApf.KNN_time,1)];
        eigtimeFull = [eigtimeFull; mean(detail.LApf.eig_time,1)];
        steptimeFull = [steptimeFull;squeeze(mean(runtimeFull,2))'];
        
        % Modify as needed depending on actual values
        xticklabel = [xticklabel, num2str(nbEig_vector(j))];
        xtick = [xtick, (i-1)*numel(nbEig_vector)+j];
    end
    groupSeparator(i) = numel(nbEig_vector)*i+0.5;
end

figure();
boxplot(RMSEFull');
title('RMSE');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);


figure();
boxplot(timeFull');
title('Runtime');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,90],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);

figure();
boxplot(steptimeFull');
title('step time');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);

figure();
boxplot(KNNtimeFull');
title('KNN time');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);

figure();
boxplot(eigtimeFull');
title('Eigen time');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);