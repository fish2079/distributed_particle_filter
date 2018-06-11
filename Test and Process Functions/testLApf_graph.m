function testLApf_graph(track, graphMethod, weightedEdge, weightedEdgeStyle)
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

warning('off','all');
% Set path to helper functions
addpath('./DynamicModels/');
addpath('./HelperFunctions/');
addpath('./MeasurementModels/');
addpath('./TrackingAlgorithms/');

m_vector = 16;%[4,6,9,16,25]; %[6,9,16,25,50,100,400,700,1000]; %[4,6,9,16,25,50,100,400,700,1000];

% sim_parameters.initialization = 'perfect';

% Graph parameters
KNN_vector = 5; %[50];%[10,25,50,100,200];
Epsilon_vector = 1/8; %[1/8,1/4,1/2];

sim_parameters.N = 1000;

% Number of random trials
sim_parameters.no_trials = 40; 

switch track
    case 2
        sim_parameters.measModel = 'bearing';
    case 3
        sim_parameters.measModel = 'range';
end

% Flag for parallel run
sim_parameters.parallel = false;

% Flag for visualizing at each time step
sim_parameters.visualizeParticles = false;

% Flag for using gossip or exact aggregate
sim_parameters.gossip = false;

% sim_parameters.KNNgraph = true;
sim_parameters.graphMethod = graphMethod;

sim_parameters.weightedEdge = weightedEdge;
if (nargin==4)
    sim_parameters.weightedEdgeStyle = weightedEdgeStyle;
else
    sim_parameters.weightedEdgeStyle = 1;
end

% Tracking algorithms are
% 1. centralized bootstrap PF: BS
% 2. distributed CSS PF 
% 3. distributed LC PF
% 4. distributed Graph PF
alg_lists = {@BSpf, @CSSpf_distributed, @LCpf_distributed, @LADelaunaypf_distributed, @ClusterDelaunaypf_distributed};
sim_parameters.algorithms = alg_lists(4);

% Select the track
sim_parameters.track = track;

% sim_parameters.max_gossip_iter = 100; 
switch graphMethod
    case 'Dummy'
    for i=1:numel(m_vector)
        sim_parameters.nbEig = m_vector(i);

        % Run the simulated track with all selected tracking algorithms 
        % Each filter uses N particles   
        [results, parameters]= runSimulatedTrack(sim_parameters);

        % Store the tracking results
        filename{i} = ['Track',num2str(track),'_LApf'];
        if (weightedEdge)
            filename{i} = [filename{i}, '_Weighted',num2str(sim_parameters.weightedEdgeStyle)];
        else
            filename{i} = [filename{i}, '_Unweighted'];
        end
        filename{i} = [filename{i},'_Dummy'];
        filename{i} = [filename{i},'_m',num2str(parameters.F.LA.m)];
        filename{i} = [filename{i},'_N',num2str(parameters.F.N)];
        filename{i} = [filename{i},'_trials',num2str(parameters.no_trials)];
        filename{i} = [filename{i},'.mat'];
        save(filename{i}, 'results','parameters');
    end
    
    case 'Delaunay'
    for i=1:numel(m_vector)
        sim_parameters.nbEig = m_vector(i);

        % Run the simulated track with all selected tracking algorithms 
        % Each filter uses N particles   
        [results, parameters]= runSimulatedTrack(sim_parameters);

        % Store the tracking results
        filename{i} = ['Track',num2str(track),'_LApf'];
        if (weightedEdge)
            filename{i} = [filename{i}, '_Weighted',num2str(sim_parameters.weightedEdgeStyle)];
        else
            filename{i} = [filename{i}, '_Unweighted'];
        end
        filename{i} = [filename{i},'_DT'];
        filename{i} = [filename{i},'_m',num2str(parameters.F.LA.m)];
        filename{i} = [filename{i},'_N',num2str(parameters.F.N)];
        filename{i} = [filename{i},'_trials',num2str(parameters.no_trials)];
        filename{i} = [filename{i},'.mat'];
        save(filename{i}, 'results','parameters');
    end

    case 'KNN'
    for j=1:numel(KNN_vector)
        sim_parameters.KNN = KNN_vector(j);
        for i=1:numel(m_vector)
            sim_parameters.nbEig = m_vector(i);

            % Run the simulated track with all selected tracking algorithms 
            % Each filter uses N particles   
            [results, parameters]= runSimulatedTrack(sim_parameters);

            % Store the tracking results
            filename{i} = ['Track',num2str(track),'_LApf'];
            if (weightedEdge)
                filename{i} = [filename{i}, '_Weighted',num2str(sim_parameters.weightedEdgeStyle)];
            else
                filename{i} = [filename{i}, '_Unweighted'];
            end
            filename{i} = [filename{i}, '_KNN',num2str(parameters.KNN)];
            filename{i} = [filename{i},'_m',num2str(parameters.F.LA.m)];
            filename{i} = [filename{i},'_N',num2str(parameters.F.N)];
            filename{i} = [filename{i},'_trials',num2str(parameters.no_trials)];
            filename{i} = [filename{i},'.mat'];
            save(filename{i}, 'results','parameters');
        end
    end

    case 'Epsilon'
    for j=1:numel(Epsilon_vector)
       sim_parameters.Epsilon = Epsilon_vector(j);
       for i=1:numel(m_vector)
           sim_parameters.nbEig = m_vector(i);
        
           % Run the simulated track with all selected tracking algorithms 
           % Each filter uses N particles   
           [results, parameters]= runSimulatedTrack(sim_parameters);

           % Store the tracking results
           filename{i} = ['Track',num2str(track),'_LApf'];
           if (weightedEdge)
               filename{i} = [filename{i}, '_Weighted',num2str(sim_parameters.weightedEdgeStyle)];
           else
               filename{i} = [filename{i}, '_Unweighted'];
           end
           filename{i} = [filename{i}, '_Epsilon',num2str(round(1/parameters.Epsilon))];
           filename{i} = [filename{i},'_m',num2str(parameters.F.LA.m)];
           filename{i} = [filename{i},'_N',num2str(parameters.F.N)];
           filename{i} = [filename{i},'_trials',num2str(parameters.no_trials)];
           filename{i} = [filename{i},'.mat'];
           save(filename{i}, 'results','parameters');
       end
   end
end
