warning('off','all');
clear;clc;

filepath = '';
% Number of particles for the filter
maxDegree_vector = [1:10];

% Number of random trials
sim_parameters.no_trials = 200; 

RMSE = zeros(1, numel(maxDegree_vector));
stdev = zeros(1, numel(maxDegree_vector));
time = zeros(1, numel(maxDegree_vector));
for i=1:numel(maxDegree_vector)
    % Set number of particles
    sim_parameters.maxDegree = maxDegree_vector(i); 

    filename{1} = [filepath,'maxDegree',num2str(sim_parameters.maxDegree),'_trials',num2str(sim_parameters.no_trials),'.mat'];

    [RMSE_vector, time_vector] = processResults(filename);

    RMSE(i) = mean(RMSE_vector);
    stdev(i) = std(RMSE_vector);
    time(i) = mean(time_vector);
end

figure();
set(gcf,'color','white');
stem(0:9, RMSE, 'linewidth',6);
xlabel('Max degree of LCpf');
ylabel('RMSE');
set(gca,'fontsize',32);
xlim([0,9]);

figure();
set(gcf,'color','white');
stem(0:9, time, 'linewidth',6);
xlabel('Max degree of LCpf');
ylabel('Total runtime');
set(gca,'fontsize',32);
xlim([0,9]);
% [x,y]=meshgrid(nbEig_vector,maxDegree_vector);
% scatter(x(:),y(:), 200, RMSE(:), 'filled');
% colorbar;
% xlabel('Number of Eigenvectors');
% ylabel('Number of nearest neighbors');
% title('Average RMSE for LA filter');
% 
% figure();
% set(gcf,'color','white');
% set(gca,'fontsize',32);
% scatter(x(:),y(:), 200, time(:), 'filled');
% colorbar;
% xlabel('Number of Eigenvectors');
% ylabel('Number of nearest neighbors');
% title('Average runtime for LA filter');
