warning('off','all');
clear;clc;

filepath = 'Results_LCpf\';
% Number of particles for the filter
maxDegree_vector = [2:10];

% Number of random trials
sim_parameters.no_trials = 200; 

RMSE = zeros(1, numel(maxDegree_vector));
time = zeros(1, numel(maxDegree_vector));

RMSEFull = [];
timeFull = [];

xticklabel = {};
for i=1:numel(maxDegree_vector)
    % Set number of particles
    sim_parameters.max_degree = maxDegree_vector(i); 

    filename{1} = [filepath,'maxDegree',num2str(sim_parameters.max_degree),'_trials',num2str(sim_parameters.no_trials),'.mat'];

    [RMSE_vector, time_vector] = processResults(filename);
    RMSEFull = [RMSEFull; RMSE_vector];

    RMSE(i) = mean(RMSE_vector);
    time(i) = mean(time_vector);
    timeFull = [timeFull; time_vector];
    
    xticklabel{i} = num2str(maxDegree_vector(i)-1);
end

figure();
set(gcf,'color','white');
boxplot(RMSEFull');
xlabel('Max degree');
ylabel('RMSE');
set(gca,'fontsize',32);
set(gca,'xtick', 1:numel(maxDegree_vector));
set(gca,'xticklabel', xticklabel);

figure();
set(gcf,'color','white');
boxplot(timeFull');
xlabel('Max degree');
ylabel('Total runtime');
set(gca,'fontsize',32);
set(gca,'xtick', 1:numel(maxDegree_vector));
set(gca,'xticklabel', xticklabel);