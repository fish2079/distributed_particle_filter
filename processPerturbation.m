% Number of particles for the filter
perturbation_vector = [1, 5,10,25, 50,100,200,300]; %[1, 5,10,25];%,[50,100,200,500];

% Number of random trials
sim_parameters.no_trials = 40; 

% Select the track
sim_parameters.track = 3;

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
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([1,numel(perturbation_vector)]);
