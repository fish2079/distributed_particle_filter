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
KNN_vector = [3:10, 20, 50];
nbEig_vector = [6,10, 100, 500];%, 1000];

% Number of random trials
sim_parameters.no_trials = 100; 

RMSE = zeros(numel(KNN_vector), numel(nbEig_vector));
time = zeros(numel(KNN_vector), numel(nbEig_vector));
timeSTDEV = time;

RMSEFull = [];
timeFull = [];

KNNtimeFull = [];
eigtimeFull = [];
steptimeFull = [];
weightdifFull = [];

% xlabel parameters for the plot
xticklabel = {};
xtick = [];
groupSeparator = [];

for j=1:numel(nbEig_vector)
    sim_parameters.nbEig = nbEig_vector(j);
    for i=1:numel(KNN_vector)
        % Set number of particles
        sim_parameters.KNN = KNN_vector(i); 

        filename{i} = [filepath,'LApf_KNN',num2str(sim_parameters.KNN),'_nbEig',num2str(sim_parameters.nbEig),'_trials',num2str(sim_parameters.no_trials),'.mat'];

        [RMSE_vector, ~, time_vector, runtimeFull, detail] = extractResults(filename{i});
        
        RMSE(i,j) = nanmean(RMSE_vector);
        RMSEFull = [RMSEFull; RMSE_vector];
        
        
        timeFull = [timeFull; squeeze(sum(runtimeFull,2))'];
        time(i,j) = mean(squeeze(sum(runtimeFull,2))');
        timeSTDEV(i,j) = std(squeeze(sum(runtimeFull,2))');
        
        KNNtimeFull = [KNNtimeFull; sum(detail.LApf.KNN_time,1)];
        eigtimeFull = [eigtimeFull; sum(detail.LApf.eig_time,1)];
        steptimeFull = [steptimeFull;squeeze(sum(runtimeFull,2))'];
        weightdifFull = [weightdifFull; mean(detail.LApf.weight_dif,1)];
        
        % Modify as needed depending on actual values
        xticklabel = [xticklabel, num2str(KNN_vector(i))];
        xtick = [xtick, (j-1)*numel(KNN_vector)+i];
    end
    groupSeparator(j) = numel(KNN_vector)*j+0.5;
end
% 
% figure();
% set(gcf,'color','white');
% bar([mean(eigtimeFull,2),mean(KNNtimeFull,2), mean(steptimeFull,2)-(mean(KNNtimeFull,2)+mean(eigtimeFull,2))],'stacked');
% xlabel('K');
% ylabel('time (s)');
% legend('Eigenvalue decomposition', 'KNN graph construction', 'Miscellaneous computation');
% set(gca,'fontsize',32);
% ylim([0,65]);
% set(gca,'xticklabel', xticklabel);

figure();
set(gcf,'color','white');
boxplot(weightdifFull');
ylabel('Approximation error');
xlabel('K');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick(1:2:end));
set(gca,'xticklabel', xticklabel(1:2:end));
set(gca,'fontsize',32);
ylim([0,0.15]);
for i=1:numel(nbEig_vector)
    dim = [.1+0.25*(i-1) .5 .3 .3];
    str = ['m = ',num2str(nbEig_vector(i))];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',32);
end

figure();
set(gcf,'color','white');
boxplot(RMSEFull');
ylabel('RMSE');
xlabel('K');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
% hold on;
% plot(0:10, 2.0577*ones(1,11),'k--','linewidth',6);
set(gca,'xtick', xtick(1:2:end));
set(gca,'xticklabel', xticklabel(1:2:end));
set(gca,'fontsize',32);
ylim([0,8]);
for i=1:numel(nbEig_vector)
    dim = [.1+0.25*(i-1) .5 .3 .3];
    str = ['m = ',num2str(nbEig_vector(i))];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',32);
end

figure();
boxplot(timeFull');
ylabel('total Runtime (s)');
xlabel('K');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,90],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);

return;
figure();
boxplot(steptimeFull');
ylabel('Runtime per time step (s)');
xlabel('K');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);

figure();
boxplot(KNNtimeFull');
ylabel('KNN graph construction time (s)');
xlabel('K');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);

figure();
boxplot(eigtimeFull');
ylabel('Eigendecomposition time (s)');
xlabel('K');
hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end
set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',32);