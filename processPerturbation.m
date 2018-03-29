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

clear;clc;

% Number of particles for the filter
perturbation_vector = [1, 5,10,25, 50,100,200,500]; %[1, 5,10,25];%,[50,100,200,500];

% Number of random trials
sim_parameters.no_trials = 200; 

sim_parameters.max_gossip_iter = 100;

% Flag for parallel run
sim_parameters.parallel = true;

% Flag for visualizing at each time step
sim_parameters.visualizeParticles = false;

% Flag for using gossip or exact aggregate
sim_parameters.gossip = false;

% Select the track
sim_parameters.track = 2;

%%
sim_parameters.nbEig = 6;
%%

N = 500;
A_error = [];
w_error = [];
% Loop through each choice of particle number
for i=1:numel(perturbation_vector)
    % Set number of particles
    sim_parameters.perturbation = perturbation_vector(i)/100;
    
   % Store the tracking results
    filename{i} = ['Track',num2str(sim_parameters.track)]; 
    filename{i} = [filename{i}, '_perturbation',num2str(sim_parameters.perturbation*100)];
    filename{i} = [filename{i},'_N',num2str(N)];
    filename{i} = [filename{i},'_trials',num2str(sim_parameters.no_trials)];
    filename{i} = [filename{i},'.mat'];

    load(filename{i});
    weight_error = [];
    AER = [];
    Neff = [];
    for tr=1:parameters.no_trials
        weight_error = cat(3, weight_error, results.details{tr}{1}.weight_error);
        AER = cat(3, AER, results.details{tr}{1}.AER);
        Neff = cat(3, Neff, results.details{tr}{1}.Neff);
    end
    w_error= [w_error; mean(mean(weight_error,3),2)'];
    A_error = [A_error; mean(mean(AER,3),2)'];
    xticklabel{i} = num2str(perturbation_vector(i));
end

w_error = w_error';

color = {'m','g','k','r','b'};
figure();
set(gcf,'color','white');
hold on;
legendText = {};
for i=1:size(w_error,1)
    plot(w_error(i,:),color{i},'linewidth',6);
end
xlabel('Perturbation');
ylabel({'Weight', 'discrepancy'});
set(gca,'xtick',1:numel(perturbation_vector));
set(gca,'xticklabel',xticklabel);
legend({'CSSpf', 'LCpf measurement', 'LCpf llh','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([1,numel(perturbation_vector)]);
