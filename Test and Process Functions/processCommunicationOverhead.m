% Test script for comparing the algorithms' performance in terms of communication
% overhead on the simulated track
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

bsRMSE = 0.8291; 

filepath = 'All PF Results\Overhead_Range\';

% Number of particles for the filter
N = 500; 
gossip_vector = [1:10,20,50,100];
% overhead_factor = [6,9,9,6,6];
overhead_factor = [9,9,6,6];

% Number of random trials
sim_parameters.no_trials = 200;

% sim_parameters.max_gossip_iter = 100;

sim_parameters.N = N; 

RMSEFull = [];
steptimeFull = [];
weight_difFull = [];
% N_effFull = [];
aggregate_error_ratioFull = [];

for i=1:numel(gossip_vector)
    RMSEFull_sf = [];
    steptimeFull_sf = [];
    weight_difFull_sf = [];
    % N_effFull = [];
    aggregate_error_ratioFull_sf = [];

    % Store the tracking results
    filename = [filepath, 'Track3_'];
    filename = [filename, '_overhead',num2str(gossip_vector(i))];
    filename = [filename, '_gossip',num2str(gossip_vector(i))];
    filename = [filename,'_N',num2str(sim_parameters.N)];
    filename = [filename,'_trials',num2str(sim_parameters.no_trials)];
    filename = [filename,'.mat'];

    data = load(filename);
    [RMSEFull_alg, steptimeFull_alg, weight_dif_full_alg, ~, aggregate_error_ratio_alg, ~] = extractResults(data.results, data.parameters, false);

    RMSEFull_sf = cat(1,RMSEFull_sf, RMSEFull_alg);
    steptimeFull_sf = cat(1, steptimeFull_sf, steptimeFull_alg);
    weight_difFull_sf = cat(1, weight_difFull_sf, mean(abs(weight_dif_full_alg),4));
    aggregate_error_ratioFull_sf = cat(1, aggregate_error_ratioFull_sf, aggregate_error_ratio_alg);

    RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
    steptimeFull = cat(4, steptimeFull, steptimeFull_sf);
    weight_difFull = cat(4, weight_difFull, mean(abs(weight_difFull_sf),4));
    aggregate_error_ratioFull = cat(4, aggregate_error_ratioFull, aggregate_error_ratioFull_sf);    
    
    xticklabel{i} = num2str(gossip_vector(i));
end

RMSEFull(1,:,:,:) = [];
steptimeFull(1,:,:,:) = [];
weight_difFull(1,:,:,:) = [];
aggregate_error_ratioFull(1,:,:,:) = [];

% plot_color = {'m','g','k','r','b','c'}; %{'y','r','b','g','k','m'};\
plot_color = {'g','k','r','b','c'}; %{'y','r','b','g','k','m'};

% figure();
% set(gcf,'color','white');
% xlabel('Communication overhead per sensor per time step');
% ylabel('Weight error ratio')
% weight_plot = (squeeze(mean(mean(weight_difFull,2),3)));
% hold on;
% for i=1:size(weight_plot,1)
%     plot(gossip_vector*overhead_factor(i), weight_plot(i,:),plot_color{i}, 'linewidth',8);
% end
% % set(gca,'xtick', 1:numel(gossip_vector));
% % set(gca,'xticklabel', xticklabel);
% set(gca,'fontsize',45);
% % legend('CSSpf', 'LCpf-meas (d=20)','LCpf-llh (d=4)','LCpf-llh (d=9)','LApf (m=6)','Clusterpf (C=6)');
% legend('CSSpf', 'LCpf-llh (d=9)','LCpf-llh-GS (d=9)','LApf (m=6)','Clusterpf (C=6)');

figure();
set(gcf,'color','white');
% xlabel('Communication overhead');
xlabel('NGossip');
ylabel('RMSE')
RMSE_plot = (squeeze(mean(mean(RMSEFull,2),3)));
hold on;
for i=1:size(RMSE_plot,1)
%     plot(gossip_vector*overhead_factor(i),RMSE_plot(i,:),strcat(plot_color{i},'o-'), 'linewidth',8,'markersize',16);
    plot(1:numel(gossip_vector),RMSE_plot(i,:),strcat(plot_color{i},'o-'), 'linewidth',8,'markersize',16);
end
plot(1:numel(gossip_vector),ones(1,numel(gossip_vector))*bsRMSE,strcat(plot_color{end},'o-'), 'linewidth',8,'markersize',16);
set(gca,'xtick', 1:numel(gossip_vector));
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',45);
% legend('CSSpf', 'LCpf-meas (d=20)','LCpf-llh (d=4)','LCpf-llh (d=9)','LApf (m=6)','Clusterpf (C=6)');
% legend('CSSpf (6)', 'LCpf (9)','LCpf-GS (9)','LApf (6)','Clusterpf (6)');
legend('LCpf (9)','LCpf-GS (9)','LApf (6)','Clusterpf (6)','BSpf');

figure();
set(gcf,'color','white');
xlabel('NGossip');
ylabel('Runtime')
runtime_plot = squeeze(sum(mean(steptimeFull,3),2));
hold on;
for i=1:size(runtime_plot,1)
    plot(1:numel(gossip_vector),runtime_plot(i,:),strcat(plot_color{i},'o-'), 'linewidth',8,'markersize',16);
end
set(gca,'xtick', 1:numel(gossip_vector));
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',45);
% legend('CSSpf', 'LCpf-meas (d=20)','LCpf-llh (d=4)','LCpf-llh (d=9)','LApf (m=6)','Clusterpf (C=6)');
% legend('CSSpf (6)', 'LCpf (9)','LCpf-GS (9)','LApf (6)','Clusterpf (6)');
legend('LCpf (9)','LCpf-GS (9)','LApf (6)','Clusterpf (6)');

figure();
set(gcf,'color','white');
xlabel('Communication overhead per sensor per time step');
ylabel('Weight Error')
weight_dif_plot = (squeeze(mean(mean(weight_difFull,2),3)));
hold on;
for i=1:size(weight_dif_plot,1)
%     plot(gossip_vector*overhead_factor(i),RMSE_plot(i,:),strcat(plot_color{i},'o-'), 'linewidth',8,'markersize',16);
    plot(1:numel(gossip_vector),weight_dif_plot(i,:),strcat(plot_color{i},'o-'), 'linewidth',8,'markersize',16);
end
set(gca,'xtick', 1:numel(gossip_vector));
set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',45);
% legend('CSSpf', 'LCpf-meas (d=20)','LCpf-llh (d=4)','LCpf-llh (d=9)','LApf (m=6)','Clusterpf (C=6)');
% legend('CSSpf (6)', 'LCpf (9)','LCpf-GS (9)','LApf (6)','Clusterpf (6)');
legend('LCpf (9)','LCpf-GS (9)','LApf (6)','Clusterpf (6)');


% figure();
% set(gcf,'color','white');
% xlabel('Communication overhead per sensor per time step');
% ylabel('AER')
% RMSE_plot = (squeeze(mean(mean(aggregate_error_ratioFull,2),3)));
% hold on;
% for i=1:size(RMSE_plot,1)
%     plot(gossip_vector*overhead_factor(i),RMSE_plot(i,:),strcat(plot_color{i},'o-'), 'linewidth',8,'markersize',16);
% end
% % set(gca,'xtick', 1:numel(gossip_vector));
% % set(gca,'xticklabel', xticklabel);
% set(gca,'fontsize',45);
% % legend('CSSpf', 'LCpf-meas (d=20)','LCpf-llh (d=4)','LCpf-llh (d=9)','LApf (m=6)','Clusterpf (C=6)');
% % legend('CSSpf', 'LCpf (d=9)','LCpf-GS (d=9)','LApf (m=6)','Clusterpf (C=6)');
% legend('CSSpf', 'LCpf (d=9)','LCpf-GS (d=9)','LApf (m=6)','Clusterpf (C=6)');