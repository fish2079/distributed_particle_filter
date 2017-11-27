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
KNN_vector = [5,10,20,50];
nbEig_vector = [50:50:1000];

% Number of random trials
sim_parameters.no_trials = 100; 

RMSE = zeros(numel(KNN_vector), numel(nbEig_vector));
time = zeros(numel(KNN_vector), numel(nbEig_vector));
for j=1:numel(nbEig_vector)
    sim_parameters.nbEig = nbEig_vector(j);
    for i=1:numel(KNN_vector)
        % Set number of particles
        sim_parameters.KNN = KNN_vector(i); 

        filename{1} = [filepath,'KNN',num2str(sim_parameters.KNN),'_nbEig',num2str(sim_parameters.nbEig),'_trials',num2str(sim_parameters.no_trials),'.mat'];

        [RMSE_vector, time_vector] = processResults(filename);
        
        RMSE(i,j) = nanmean(RMSE_vector);
        time(i,j) = mean(time_vector);
        yo=5;
    end
end

figure();
set(gcf,'color','white');
[x,y]=meshgrid(nbEig_vector,KNN_vector);
scatter(x(:),y(:), 200, RMSE(:), 'filled');
colorbar;
xlabel('Number of Eigenvectors');
ylabel('Number of neighbors');
title('Average RMSE for LA filter');
set(gca,'fontsize',32);

figure();
set(gcf,'color','white');
scatter(x(:),y(:), 200, time(:), 'filled');
colorbar;
xlabel('Number of Eigenvectors');
ylabel('Number of neighbors');
title('Average runtime for LA filter');
set(gca,'fontsize',32);