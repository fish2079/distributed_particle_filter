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

filepath = 'All PF Results\';

% Number of particles for the filter
N = 500; 
overhead_vector = [60:60:1200];

% Number of random trials
sim_parameters.no_trials = 200;

% sim_parameters.max_gossip_iter = 100;

% Flag for parallel run
sim_parameters.parallel = true;

% Flag for visualizing at each time step
sim_parameters.visualizeParticles = false;

% Flag for using gossip or exact aggregate
sim_parameters.gossip = true;

% Tracking algorithms are
% 1. centralized bootstrap PF: BS
% 2. distributed CSS PF 
% 3. distributed LC PF
% 4. distributed Graph PF
% alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LApf_distributed, @LADelaunaypf_distributed, @Clusterpf_distributed, @ClusterDelaunaypf_distributed};
alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LADelaunaypf_distributed, @ClusterDelaunaypf_distributed};
overhead_factor = [1,6,20,6,6];
% sim_parameters.algorithms = alg_lists([1:5]);

% sim_parameters.areaLength = 125;
% Loop through each choice of particle number
sim_parameters.nbEig = 6;
sim_parameters.nbClusters = 6;

sim_parameters.N = N; 

RMSEFull = [];
steptimeFull = [];
weight_difFull = [];
% N_effFull = [];
aggregate_error_ratioFull = [];

for i=1:numel(overhead_vector)
    RMSEFull_sf = [];
    steptimeFull_sf = [];
    weight_difFull_sf = [];
    % N_effFull = [];
    aggregate_error_ratioFull_sf = [];
    for alg=2:5
        sim_parameters.algorithms = alg_lists(alg);
        sim_parameters.max_gossip_iter = overhead_vector(i)/overhead_factor(alg);

        % Store the tracking results
        filename = [filepath, 'Track3_'];
        filename = [filename, func2str(alg_lists{alg})];
        filename = [filename, '_overhead',num2str(overhead_vector(i))];
        filename = [filename, '_gossip',num2str(sim_parameters.max_gossip_iter)];
        filename = [filename,'_N',num2str(sim_parameters.N)];
        filename = [filename,'_trials',num2str(sim_parameters.no_trials)];
        filename = [filename,'.mat'];
        
        data = load(filename);
        [RMSEFull_alg, steptimeFull_alg, weight_dif_full_alg, ~, aggregate_error_ratio_alg, ~] = extractResults(data.results, data.parameters, false);

        RMSEFull_sf = cat(1,RMSEFull_sf, RMSEFull_alg);
        steptimeFull_sf = cat(1, steptimeFull_sf, steptimeFull_alg);
        weight_difFull_sf = cat(1, weight_difFull_sf, mean(abs(weight_dif_full_alg),4));
        aggregate_error_ratioFull_sf = cat(1, aggregate_error_ratioFull_sf, aggregate_error_ratio_alg);
    end
    
    RMSEFull = cat(4,RMSEFull, RMSEFull_sf);
    steptimeFull = cat(4, steptimeFull, steptimeFull_sf);
    weight_difFull = cat(4, weight_difFull, mean(abs(weight_difFull_sf),4));
    aggregate_error_ratioFull = cat(4, aggregate_error_ratioFull, aggregate_error_ratioFull_sf);    
end