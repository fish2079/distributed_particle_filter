filepath = 'All PF Results/Perturbation_Bearing/';

% Number of particles for the filter
gossip_vector = [1:10,20,30, 50,75, 100, 125];

% Number of random trials
sim_parameters.no_trials = 40; 

% Select the track
sim_parameters.track = 2;

N = 500;
A_error = [];
w_error = [];
n_error = [];
n_error_time = [];
Neff_full = [];

% Loop through each choice of particle number
for i=1:numel(gossip_vector)
    % Set number of gossip iterations
    sim_parameters.max_gossip_iter = gossip_vector(i);
    
   % Store the tracking results
    filename{i} = [filepath,'Track',num2str(sim_parameters.track)]; 
    filename{i} = [filename{i}, '_Gossip',num2str(sim_parameters.max_gossip_iter)];
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
%         results.details{tr}{1}.errorNorm(:,7,:) = max(results.details{tr}{1}.errorNorm(:,1,:),results.details{tr}{1}.errorNorm(:,2,:));
        results.details{tr}{1}.errorNorm(:,13,:) = results.details{tr}{1}.errorNorm(:,11,:).*results.details{tr}{1}.errorNorm(:,12,:);
        normError = cat(4, normError, results.details{tr}{1}.errorNorm);
    end
    w_error= [w_error; trimmean(trimmean(weight_error,5,3),5,2)'];
    A_error = [A_error; trimmean(trimmean(AER,5,3),5,2)'];
    Neff_full = cat(3, Neff_full, mean(Neff,3));
    n_error = cat(3, n_error, trimmean(trimmean(normError,5,4),5,3));
    n_error_time = cat(4, n_error_time, trimmean(normError,5,4));
    xticklabel{i} = num2str(gossip_vector(i));
end

w_error = w_error';

color = {'m','g','k','r','b'};
% color = {'g','k','r','b'};

% n_error_time(1,:,:,:) = [];
% n_error(1,:,:) = [];
% w_error(1,:) = [];

% encoding error over time
delta_m_time = squeeze(n_error_time(:,1,:,2));
figure();
set(gcf,'color','white');
hold on;
for i=1:size(delta_m_time,1)
    plot(delta_m_time(i,:),color{i},'linewidth',6);
end
xlabel('Time');
ylabel({'\delta_m'});
set(gca,'YScale','log');
plot(ones(1,numel(delta_m_time(i,:))),'c:','linewidth',6);
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% legend({'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([2,50]);

% gossip error wrt NGossip
delta_gossip = squeeze(n_error(:,2,:));
delta_gossip_limit = squeeze(n_error(:,8,:));
delta_gossip_limit2 = squeeze(n_error(:,9,:));
figure();
set(gcf,'color','white');
hold on;
for i=1:size(delta_gossip,1)
    plot(delta_gossip(i,:),strcat('o',color{i}),'linewidth',6,'markersize',16);
end
for i=1:size(delta_gossip,1)
    plot(delta_gossip_limit(i,:),color{i},'linewidth',6);
end
% for i=1:size(delta_gossip,1)
%     plot(delta_gossip_limit2(i,:),strcat('--',color{i}),'linewidth',6);
% end
xlabel('NGossip');
ylabel({'\delta_{gossip}'});
set(gca,'xtick',1:numel(gossip_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
plot(ones(1,numel(gossip_vector)),'c:','linewidth',6);
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% legend({'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([1,numel(gossip_vector)]);

delta_gossip = squeeze(n_error(:,2,:));
delta_gossip_limit = squeeze(n_error(:,8,:));
delta_gossip_limit2 = squeeze(n_error(:,9,:));
figure();
set(gcf,'color','white');
hold on;
for i=1:size(delta_gossip,1)
    plot(delta_gossip(i,:),strcat('o',color{i}),'linewidth',6,'markersize',16);
end
% for i=1:size(delta_gossip,1)
%     plot(delta_gossip_limit(i,:),color{i},'linewidth',6);
% end
for i=1:size(delta_gossip,1)
    plot(delta_gossip_limit2(i,:),strcat('--',color{i}),'linewidth',6);
end
xlabel('NGossip');
ylabel({'\delta_{gossip}'});
set(gca,'xtick',1:numel(gossip_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
plot(ones(1,numel(gossip_vector)),'c:','linewidth',6);
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% legend({'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([1,numel(gossip_vector)]);

% Beta plot
Beta = squeeze(n_error(:,3,:));
figure();
set(gcf,'color','white');
hold on;
for i=1:size(Beta,1)
    plot(Beta(i,:),color{i},'linewidth',6);
end
xlabel('NGossip');
ylabel('\beta');
set(gca,'xtick',1:numel(gossip_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
plot(ones(1,numel(gossip_vector)),'c:','linewidth',6);
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% legend({'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([1,numel(gossip_vector)]);


% Plot for delta
delta = squeeze(n_error(:,5,:));
delta_limit_upper = squeeze(n_error(:,6,:));
figure();
set(gcf,'color','white');
hold on;
for i=1:size(delta,1)
    plot(delta(i,:),strcat('o',color{i}),'linewidth',6,'markersize',16);
end
for i=1:size(delta,1)
    plot(delta_limit_upper(i,:),color{i},'linewidth',6);
end
xlabel('NGossip');
ylabel({'\delta'});
set(gca,'xtick',1:numel(gossip_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
plot(ones(1,numel(gossip_vector)),'c:','linewidth',6);
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% legend({'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([1,numel(gossip_vector)]);

% Plot for exp(\Psi \delta \alpha)$
% Baha = squeeze(n_error(:,7,:));
% figure();
% set(gcf,'color','white');
% hold on;
% for i=1:size(Baha,1)
%     plot(Baha(i,:),color{i},'linewidth',6);
% end
% xlabel('NGossip');
% ylabel('$$\min_i\exp(\Psi_i\hat{\alpha}-\Psi_i\alpha)$$','Interpreter','latex');
% set(gca,'xtick',1:numel(gossip_vector));
% set(gca,'xticklabel',xticklabel);
% set(gca,'YScale','log');
% plot(ones(1,numel(gossip_vector)),'c:','linewidth',6);
% legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% set(gca,'fontsize',45);
% xlim([1,numel(gossip_vector)]);

% Plot for \Psi
% Psi = squeeze(n_error_time(:,11,:,1).*n_error_time(:,12,1));
Psi = squeeze(n_error_time(:,13,:,2));
figure();
set(gcf,'color','white');
hold on;
for i=1:size(Psi,1)
    plot(Psi(i,:),color{i},'linewidth',6);
end
xlabel('Time');
ylabel({'||\Psi||_{\infty} ||\alpha||_{\infty}'});
set(gca,'xtick',1:size(Psi,2));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf'});
% legend({'LCpf', 'LCpf-GS','LApf'});
set(gca,'fontsize',45);

% Plot for \alpha
% alpha = squeeze(n_error_time(:,12,:,1));
% figure();
% set(gcf,'color','white');
% hold on;
% for i=1:size(alpha,1)
%     plot(alpha(i,:),color{i},'linewidth',6);
% end
% xlabel('Time');
% ylabel({'Alpha'});
% set(gca,'xtick',1:size(alpha,2));
% set(gca,'xticklabel',xticklabel);
% set(gca,'YScale','log');
% plot(ones(1,50),'c:','linewidth',6);
% legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% set(gca,'fontsize',45);


figure();
set(gcf,'color','white');
hold on;
legendText = {};
for i=1:size(w_error,1)
    plot(w_error(i,:),color{i},'linewidth',6);
end
xlabel('NGossip');
ylabel({'Weight', 'discrepancy'});
set(gca,'xtick',1:numel(gossip_vector));
set(gca,'xticklabel',xticklabel);
% set(gca,'YScale','log');
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% legend({'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([1,numel(gossip_vector)]);

% figure();
% set(gcf,'color','white');
% hold on;
% for i=1:size(w_error,1)
%     plot(squeeze(n_error(i,4,:))',color{i},'linewidth',6);
% end
% legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% xlabel('NGossip');
% ylabel({'1/N |1-\exp(\Psi\Delta)|'});
% set(gca,'xtick',1:numel(gossip_vector));
% set(gca,'xticklabel',xticklabel);
% set(gca,'YScale','log');
% set(gca,'fontsize',45);
% xlim([1,numel(gossip_vector)]);

% figure();
% set(gcf,'color','white');
% hold on;
% for i=1:size(w_error,1)
%     plot(squeeze(n_error(i,7,:))',color{i},'linewidth',6);
% end
% legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% xlabel('NGossip');
% ylabel({'||\Psi\Delta||_2'});
% set(gca,'xtick',1:numel(gossip_vector));
% set(gca,'xticklabel',xticklabel);
% set(gca,'fontsize',45);
% set(gca,'YScale','log');
% xlim([1,numel(gossip_vector)]);