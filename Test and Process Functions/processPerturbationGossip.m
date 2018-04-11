% Number of particles for the filter
gossip_vector = [1:10,20,50,100];

% Number of random trials
sim_parameters.no_trials = 40; 

% Select the track
sim_parameters.track = 2;

N = 500;
A_error = [];
w_error = [];
n_error = [];
n_error_time = [];

% Loop through each choice of particle number
for i=1:numel(gossip_vector)
    % Set number of gossip iterations
    sim_parameters.max_gossip_iter = gossip_vector(i);
    
   % Store the tracking results
    filename{i} = ['All PF Results/Track',num2str(sim_parameters.track)]; 
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
        results.details{tr}{1}.errorNorm(:,7,:) = max(results.details{tr}{1}.errorNorm(:,1,:),results.details{tr}{1}.errorNorm(:,2,:));
        normError = cat(4, normError, results.details{tr}{1}.errorNorm(:,[1,2,4,5,6,7],:));
    end
    w_error= [w_error; mean(mean(weight_error,3),2)'];
    A_error = [A_error; mean(mean(AER,3),2)'];
    n_error = cat(3, n_error, trimmean(mean(normError,4),5,3));
    n_error_time = cat(4, n_error_time,mean(normError,4));
    xticklabel{i} = num2str(gossip_vector(i));
end

w_error = w_error';

color = {'m','g','k','r','b'};

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
set(gca,'fontsize',45);
xlim([2,50]);

delta_gossip_time_1 = squeeze(n_error_time(:,2,:,3));
delta_gossip_limit_time_1 = squeeze(n_error_time(:,3,:,3));
figure();
set(gcf,'color','white');
hold on;
for i=1:size(delta_gossip_time_1,1)
    plot(delta_gossip_time_1(i,:),color{i},'linewidth',6);
end
for i=1:size(delta_gossip_time_1,1)
    plot(delta_gossip_limit_time_1(i,:),strcat('-o',color{i}),'linewidth',6);
end
xlabel('Time');
ylabel({'\delta_{gossip}'});
set(gca,'YScale','log');
plot(ones(1,numel(delta_gossip_time_1(i,:))),'k--','linewidth',6);
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);

% delta_gossip_time_high = squeeze(n_error_time(:,2,:,11));
% delta_gossip_limit_time_high = squeeze(n_error_time(:,3,:,11));
% figure();
% set(gcf,'color','white');
% hold on;
% for i=1:size(delta_gossip_time_high,1)
%     plot(delta_gossip_time_high(i,:),color{i},'linewidth',6);
% end
% for i=1:size(delta_gossip_time_high,1)
%     plot(delta_gossip_limit_time_high(i,:),strcat('-o',color{i}),'linewidth',6);
% end
% xlabel('Time');
% ylabel({'\delta_{gossip}'});
% set(gca,'YScale','log');
% plot(ones(1,numel(delta_gossip_time_high(i,:))),'k--','linewidth',6);
% legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% set(gca,'fontsize',45);

delta_gossip = squeeze(n_error(:,2,:));
delta_gossip_limit = squeeze(n_error(:,3,:));
figure();
set(gcf,'color','white');
hold on;
for i=1:size(delta_gossip,1)
    plot(delta_gossip(i,:),color{i},'linewidth',6);
end
for i=1:size(delta_gossip,1)
    plot(delta_gossip_limit(i,:),strcat('-o',color{i}),'linewidth',6,'markersize',32);
end
xlabel('NGossip');
ylabel({'\delta_{gossip}'});
set(gca,'xtick',1:numel(gossip_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
plot(ones(1,numel(gossip_vector)),'c:','linewidth',6);
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([1,numel(gossip_vector)]);

delta = squeeze(n_error(:,4,:));
delta_limit_upper = squeeze(n_error(:,5,:));
delta_limit_lower = squeeze(n_error(:,6,:));
figure();
set(gcf,'color','white');
hold on;
for i=1:size(delta,1)
    plot(delta(i,:),color{i},'linewidth',6);
end
for i=1:size(delta,1)
    plot(delta_limit_upper(i,:),strcat('--o',color{i}),'linewidth',6,'markersize',32);
end
xlabel('NGossip');
ylabel({'\delta'});
set(gca,'xtick',1:numel(gossip_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
plot(ones(1,numel(gossip_vector)),'c:','linewidth',6);
legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
set(gca,'fontsize',45);
xlim([1,numel(gossip_vector)]);

% figure();
% set(gcf,'color','white');
% hold on;
% legendText = {};
% for i=1:size(w_error,1)
%     plot(w_error(i,:),color{i},'linewidth',6);
% end
% xlabel('NGossip');
% ylabel({'Weight', 'discrepancy'});
% set(gca,'xtick',1:numel(gossip_vector));
% set(gca,'xticklabel',xticklabel);
% % set(gca,'YScale','log');
% legend({'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'});
% set(gca,'fontsize',45);
% xlim([1,numel(gossip_vector)]);

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