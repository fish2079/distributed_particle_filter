% Number of particles for the filter
perturbation_vector = [1, 5,10,25, 50,100,200,300]; %[1, 5,10,25];%,[50,100,200,500];

% Number of random trials
sim_parameters.no_trials = 40; 

% Select the track
sim_parameters.track = 3;

N = 500;
A_error = [];
w_error = [];
n_error = [];
Psi_time = [];

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
    normError = [];
    for tr=1:parameters.no_trials
        weight_error = cat(3, weight_error, results.details{tr}{1}.weight_error);
        AER = cat(3, AER, results.details{tr}{1}.AER);
        Neff = cat(3, Neff, results.details{tr}{1}.Neff);
        normError = cat(4, normError, results.details{tr}{1}.errorNorm);
    end
    yo = mean(normError,4);
    Psi_time = cat(4, Psi_time, squeeze(yo(:,5,:)));
    w_error= [w_error; mean(mean(weight_error,3),2)'];
    A_error = [A_error; mean(mean(AER,3),2)'];
    n_error = cat(3, n_error, mean(mean(normError,4),3));
    xticklabel{i} = num2str(perturbation_vector(i));
end

Psi_time = squeeze(Psi_time);
Psi_time = Psi_time(:,:,1);

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

% figure();
% set(gcf,'color','white');
% hold on;
% for i=1:2 %1:size(w_error,1)
%     plot(squeeze(n_error(i,1,:))',color{i},'linewidth',6);
% end
% legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% xlabel('Perturbation');
% ylabel({'||\Psi||'});
% set(gca,'xtick',1:numel(perturbation_vector));
% set(gca,'xticklabel',xticklabel);
% set(gca,'fontsize',45);
% set(gca,'YScale','log');
% xlim([1,numel(perturbation_vector)]);

figure();
set(gcf,'color','white');
hold on;
for i=1:size(w_error,1)
    plot(squeeze(n_error(i,4,:))',color{i},'linewidth',6);
end
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
xlabel('Perturbation');
ylabel({'1/N |1-exp(\Psi\Delta|)'});
set(gca,'xtick',1:numel(perturbation_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'fontsize',45);
set(gca,'YScale','log');
xlim([1,numel(perturbation_vector)]);

figure();
set(gcf,'color','white');
hold on;
for i=1:size(w_error,1)
    plot(squeeze(n_error(i,7,:))',color{i},'linewidth',6);
end
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
xlabel('Perturbation');
ylabel({'||\Psi\Delta||_2'});
set(gca,'xtick',1:numel(perturbation_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'fontsize',45);
set(gca,'YScale','log');
xlim([1,numel(perturbation_vector)]);


% figure();
% set(gcf,'color','white');
% hold on;
% for i=1:size(w_error,1)
%     plot(squeeze(n_error(i,2,:))',color{i},'linewidth',6);
% end
% legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% xlabel('Perturbation');
% ylabel({'||\Delta||_2'});
% set(gca,'xtick',1:numel(perturbation_vector));
% set(gca,'xticklabel',xticklabel);
% set(gca,'fontsize',45);
% xlim([1,numel(perturbation_vector)]);


% figure();
% set(gcf,'color','white');
% hold on;
% for i=[1,2,5]
%     plot(squeeze(Psi_time(i,:))',color{i},'linewidth',6);
% end
% legend({'CSSpf', 'LCpf','Clusterpf'});
% xlabel('time');
% ylabel({'||\Psi (\Psi^T\Psi)^{-1}||_2'});
% % set(gca,'xtick',1:numel(perturbation_vector));
% % set(gca,'xticklabel',xticklabel);
% set(gca,'YScale','log');
% set(gca,'fontsize',45);