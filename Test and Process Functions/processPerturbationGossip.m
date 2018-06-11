function processPerturbationGossip(track)
filepath = 'New Initialization Full Simulation Results\All_PF_Error\';

% Number of particles for the filter
gossip_vector = [1,5,10,25,50,75,100,200];

% Number of random trials
sim_parameters.no_trials = 200; 

N = 1000;
% A_error = [];
w_error = [];
n_error = [];
n_error_time = [];
% Neff_full = [];

% Loop through each choice of particle number
for i=1:numel(gossip_vector)    
    % Store the tracking results
    filename{i} = [filepath,'Track',num2str(track),'_Error']; 
    filename{i} = [filename{i}, '_Gossip',num2str(gossip_vector(i))];
    filename{i} = [filename{i},'_N',num2str(N)];
    filename{i} = [filename{i},'_trials',num2str(sim_parameters.no_trials)];
    filename{i} = [filename{i},'.mat'];

    load(filename{i});
    weight_error = [];
%     AER = [];
%     Neff = [];
    normError = [];
    for tr=1:parameters.no_trials
        weight_error = cat(3, weight_error, results.details{tr}{1}.weight_error);
%         AER = cat(3, AER, results.details{tr}{1}.AER);
%         Neff = cat(3, Neff, results.details{tr}{1}.Neff);
        results.details{tr}{1}.errorNorm(:,13,:) = results.details{tr}{1}.errorNorm(:,11,:).*results.details{tr}{1}.errorNorm(:,12,:);
        normError = cat(4, normError, results.details{tr}{1}.errorNorm);
    end
%     w_error= [w_error; trimmean(trimmean(weight_error,5,3),5,2)'];
    w_error= [w_error; mean(mean(weight_error,3),2)'];
%     A_error = [A_error; trimmean(trimmean(AER,5,3),5,2)'];
%     Neff_full = cat(3, Neff_full, mean(Neff,3));
%     n_error = cat(3, n_error, trimmean(trimmean(normError,5,4),5,3));
    n_error = cat(3, n_error, mean(mean(normError,4),3));
%     n_error_time = cat(4, n_error_time, trimmean(normError,5,4));
    n_error_time = cat(4, n_error_time, mean(normError,4));
    xticklabel{i} = num2str(gossip_vector(i));
end
w_error = w_error';

if (track==2)
    legendText = {'CSSpf', 'LCpf', 'LCpf-GS','LApf','Clusterpf'};
    start_index = 1;
else
    legendText = {'LCpf', 'LCpf-GS','LApf','Clusterpf'};
    start_index = 2;
end

% Figure 1
% encoding error over time
delta_m_time = squeeze(n_error_time(:,1,:,2));
figure();
set(gcf,'color','white');
hold on;
for i=start_index:size(delta_m_time,1)
    plot(delta_m_time(i,:),'linewidth',6);
end
xlabel('Time');
ylabel({'\delta_m'});
set(gca,'YScale','log');
plot(ones(1,numel(delta_m_time(i,:))),'c:','linewidth',6);
legend(legendText);
set(gca,'fontsize',32);
xlim([1,50]);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

% Fig 2 delta gossip
delta_gossip = squeeze(n_error(:,2,:));
% delta_gossip_limit = squeeze(n_error(:,8,:));
delta_gossip_limit2 = squeeze(n_error(:,9,:));
figure();
set(gcf,'color','white');
hold on;
for i=start_index:size(delta_gossip,1)
    plot(delta_gossip(i,:),'o','linewidth',6,'markersize',16);
end
set(gca,'ColorOrderIndex',1);
for i=start_index:size(delta_gossip,1)
    plot(delta_gossip_limit2(i,:),'--','linewidth',6);
end
xlabel('NGossip');
ylabel({'\delta_{gossip}'});
set(gca,'xtick',1:numel(gossip_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
plot(ones(1,numel(gossip_vector)),'c:','linewidth',6);
legend(legendText);
set(gca,'fontsize',32);
xlim([1,numel(gossip_vector)]);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

% Beta plot
Beta = squeeze(n_error(:,3,:));
figure();
set(gcf,'color','white');
hold on;
for i=start_index:size(Beta,1)-1
    plot(Beta(i,:),'linewidth',6);
end
xlabel('NGossip');
ylabel('\beta');
set(gca,'xtick',1:numel(gossip_vector));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
plot(ones(1,numel(gossip_vector)),'c:','linewidth',6);
legend(legendText(1:end-1));
set(gca,'fontsize',32);
xlim([1,numel(gossip_vector)]);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

% Plot for delta
overhead_vector = [6,4,4,9,9];
delta = squeeze(n_error(:,5,:));
delta_limit_upper = squeeze(n_error(:,6,:));
figure();
set(gcf,'color','white');
hold on;
for i=start_index:size(delta,1)
    plot(gossip_vector*overhead_vector(i), delta(i,:),'o','linewidth',6,'markersize',16);
end
set(gca,'ColorOrderIndex',1);
for i=start_index:size(delta,1)
    plot(gossip_vector*overhead_vector(i), delta_limit_upper(i,:),'linewidth',6);
end
xlabel('Scalars transmitted per sensor');
ylabel({'\delta'});
set(gca,'YScale','log');
plot(ones(1,500),'c:','linewidth',6);
legend(legendText);
set(gca,'fontsize',32);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

% Plot for \Psi
Psi = squeeze(n_error_time(:,11,:,2));
figure();
set(gcf,'color','white');
hold on;
for i=start_index:size(Psi,1)
    plot(Psi(i,:),'linewidth',6);
end
xlabel('Time');
ylabel({'||\Psi||_{\infty}'});
set(gca,'xtick',1:size(Psi,2));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
legend(legendText(1:end-1));
set(gca,'fontsize',32);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

%Plot for \alpha
Psi = squeeze(n_error_time(:,12,:,2));
figure();
set(gcf,'color','white');
hold on;
for i=start_index:size(Psi,1)
    plot(Psi(i,:),'linewidth',6);
end
xlabel('Time');
ylabel({'||\alpha||_{\infty}'});
set(gca,'xtick',1:size(Psi,2));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
legend(legendText(1:end-1));
set(gca,'fontsize',32);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

%Plot for \Psi\alpha
Psi = squeeze(n_error_time(:,13,:,2));
figure();
set(gcf,'color','white');
hold on;
for i=start_index:size(Psi,1)
    plot(Psi(i,:),'linewidth',6);
end
xlabel('Time');
ylabel({'||\Psi||_{\infty} ||\alpha||_{\infty}'});
set(gca,'xtick',1:size(Psi,2));
set(gca,'xticklabel',xticklabel);
set(gca,'YScale','log');
legend(legendText(1:end-1));
set(gca,'fontsize',32);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')

% Weight discrepancy
figure();
set(gcf,'color','white');
hold on;
for i=start_index:size(w_error,1)
    plot(gossip_vector*overhead_vector(i), w_error(i,:),'linewidth',6);
end
xlabel('NGossip');
ylabel({'Weight', 'discrepancy'});
set(gca,'YScale','log');
legend(legendText);
set(gca,'fontsize',32);
grid on;
set(gca,'GridAlpha',0.5)
set(gca,'GridColor','k')