function processCommunicationOverhead(track)
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

bsRMSE = 0.9369;
%bsRMSE = 1.6820; %1.5587; %0.8238; 

filepath = 'New Initialization Full Simulation Results\All_PF_RMSE\';

% Number of particles for the filter
N = 1000; 
gossip_vector = [1,5,10,25,50,75,100,200];
overhead_factor = [6,9,9,9,25,14];

% Number of random trials
sim_parameters.no_trials = 200;

sim_parameters.N = N; 

RMSEFull = [];
% steptimeFull = [];
% weight_difFull = [];
% N_effFull = [];
% aggregate_error_ratioFull = [];

for i=1:numel(gossip_vector)
    RMSEFull_sf = [];

    % Store the tracking results
    filename = [filepath, 'Track',num2str(track),'_Allpf_'];
    filename = [filename, '_gossip',num2str(gossip_vector(i))];
    filename = [filename,'_N',num2str(sim_parameters.N)];
    filename = [filename,'_trials',num2str(sim_parameters.no_trials)];
    filename = [filename,'.mat'];

    data = load(filename);
    [RMSEFull_alg] = extractResults(data.results, data.parameters, false);

    RMSEFull_sf = cat(1,RMSEFull_sf, RMSEFull_alg);
    
    RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
    
%     xticklabel{i} = num2str(gossip_vector(i));
end

figure();
set(gcf,'color','white');
xlabel('Scalars transmitted per sensor');
ylabel('RMSE')
RMSE_plot = (squeeze(mean(mean(RMSEFull,2),3)));
hold on;
if (track==2)
    set(gca,'ColorOrderIndex',1);
    for i=1:size(RMSE_plot,1)
        plot(gossip_vector*overhead_factor(i),RMSE_plot(i,:),'linewidth',8);
    %     plot(1:numel(gossip_vector),RMSE_plot(i,:),strcat(plot_color{i},'o-'), 'linewidth',8,'markersize',16);
    end
else
    set(gca,'ColorOrderIndex',2);
    for i=2:size(RMSE_plot,1)
        plot(gossip_vector*overhead_factor(i),RMSE_plot(i,:),'linewidth',8);
    %     plot(1:numel(gossip_vector),RMSE_plot(i,:),strcat(plot_color{i},'o-'), 'linewidth',8,'markersize',16);
    end
end
% plot(1:2000,ones(1,2000)*bsRMSE,strcat('m','--'), 'linewidth',8,'markersize',16);
set(gca,'YScale','log');
% set(gca,'xtick', 1:numel(gossip_vector));
% set(gca,'xticklabel', xticklabel);
set(gca,'fontsize',45);
if (track==2)
    legend('CSSpf (6)', 'LCpf (9)','LCpf-GS (9)','LApf (9)','Clusterpf (25)','GApf (14)','BSpf');
else
    legend('LCpf (9)','LCpf-GS (9)','LApf (9)','Clusterpf (25)','GApf (14)','BSpf');
end
xlim([0,500]);
grid on;
set(gca,'GridAlpha',1)
set(gca,'GridColor','k')