warning('off','all');
clear;clc;

filepath = 'Results_LCpf\';
% filepath = '';
% Number of particles for the filter
maxDegree_vector = [1:9];

% Number of random trials
sim_parameters.no_trials = 200; 

RMSE = zeros(1, numel(maxDegree_vector));
RMSE_full = [];
stdev = zeros(1, numel(maxDegree_vector));
time = zeros(1, numel(maxDegree_vector));
time_full = [];
Hx_ss_dif_full = [];
for i=1:numel(maxDegree_vector)
    % Set number of particles
    sim_parameters.maxDegree = maxDegree_vector(i); 

    filename{i} = [filepath,'LCpf_maxDegree',num2str(sim_parameters.maxDegree),'_trials',num2str(sim_parameters.no_trials),'.mat'];

    [RMSE_vector, ~, time_vector, runtimeFull, detail] = extractResults(filename{i});
    
    RMSE(i) = mean(RMSE_vector);
    RMSE_full = [RMSE_full; RMSE_vector];
    stdev(i) = std(RMSE_vector);
    time(i) = mean(time_vector);
    time_full = [time_full; time_vector];
    
    Hx_ss_dif_full = [Hx_ss_dif_full; mean(sum(detail.LCpf.Hx_ss_dif,3),1)];
end

figure();
set(gcf,'color','white');
boxplot(RMSE_full');
xlabel('Max degree of LCpf');
ylabel('RMSE');
set(gca,'fontsize',32);
hold on;
plot(0:10, 1.4980*ones(1,11),'r','linewidth',6);

figure();
set(gcf,'color','white');
boxplot(time_full');
xlabel('Max degree of LCpf');
ylabel('Total runtime');
set(gca,'fontsize',32);

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
