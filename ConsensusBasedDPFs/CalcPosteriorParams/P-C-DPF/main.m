function main

clc;
clear all;
disp(fix(clock));
system('uname -n');
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock))); 
 
global parset; %global structure containing all simulation parameters

%generated offline by generateTopology and loaded:
load 'adjMtx';
load 'sensorsXY';
parset.sensorsPos = sensorsXY; %physical positions of the sensors  
parset.adjacencyMat = sparse(adjMtx);  %adjacency matrix of the graph representation of the sensor network (with communication links)
clear adjMtx;
clear sensorsXY;

if (exist('systematicResamp','file')~=3) %don't compile if the file exists
    mex systematicResamp.c;
end

%Read parset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen('parset.txt','r');
parameterRead = 1; %the ID / position of the parameter that is being read from parset

while 1
    
    line = fgetl(fid);
    
    if ~isempty(line) %skip empty lines
        
            if ~ischar(line) %if end of file is reached, tline will contain -1 -> ischar will be false and while loop terminates
                break;
            end
        
        if line(1) ~='#' %skip lines starting with #
            
            value = str2double(line); %convert the string to a number 
            
            %the switch statement must correspond to the ordering of
            %parameters in parset
            switch parameterRead
                case 1
                    parset.simTime = value; %number of discrete time instants (over which state x_n evolves)
                case 2
                    parset.avgCycles = value; %number of cycles over which we calculate (average) the MSE                    
                case 3
                    parset.simAreaSize = value; %the simulation area is a square of simAreaSize x simAreaSize [m]; this value must be the same as the simAreaSize value found in generateTopology!!!                
                case 4
                    %parset.Cu = diag([value,value,value,value]); %covariance of the Gaussian driving noise
                    parset.Cu = eye(2)*value; %covariance of the Gaussian driving noise
                case 5
                    parset.numTargets = value; %number of targets
                case 6
                    parset.A = value; %amplitude of the sound source (target)
                case 7
                    parset.w_sigma = sqrt(value); %standard deviation of the measurement noise
                case 8
                    parset.measModel = value; %measurement model
                    %0 - Inverse decay h(x) = a/(r^invPow)
                    %1 - Exponential decay h(x) = a*exp(-lambda*r)
                    %2 - Distance squared h(x) = r^2
                case 9
                    parset.invPow = value; %inverse decay parameter
                case 10
                    parset.lambda = value; %exponential decay parameter                  
                case 11
                    parset.numParticles = value; %number of particles                 
                case 12
                    parset.psProposal = value; %prediction step proposal density: 0 - state transition; 1 - auxiliary particle filtering
                case 13
                    parset.consensusType = value;%consensus type
                    %0 - perfect consensus (consensus sums are calculated exactly; asymptotic case)
                    %1 - realistic consensus
                case 14
                    if (parset.consensusType == 1)
                        parset.consensusAlgorithm = value;
                    end
                case 15
                    %expr = '(?<con_w>.*?)\[(?<b>.*?)\]';
                    %expr = '(\d+\.\d+)|(\d+)';
                    expr = '(?<type>\d+?)\[(?<opt>(\d+)?)\]|(?<type>\d+?)\[(?<opt>\d+\.\d+?)\]|(?<type>(\d+)?)';
                    [dummy,tmp2,dummy] = regexp(line, expr,'match', 'names'); %#ok<NASGU,ASGLU>
                    
                    if (isempty(tmp2.opt)==1)
                        parset.consensusWeightModel = str2double(tmp2.type);
                    else
                        parset.consensusWeightModel = [str2double(tmp2.type);str2double(tmp2.opt)];
                    end
                case 16
                    parset.consensusNumOfIter = value;
                case 17
                    %Network topology generation
                    parset.dynGenerateTopology = value; %0 - load pre-generated network topology; 1 - randomly generate topology in each avg cycle
                case 18
                    %Target trajectory generation
                    parset.pregenTrajectory = value; %0 - random generation during sim run; 1 - pre-generated target trajectory
                otherwise
                    break;       
            end
        
            parameterRead = parameterRead + 1;
        end
    end
end

fclose(fid);

[numSensors] = size(parset.sensorsPos,1); %determine the number of sensors
parset.numSensors = numSensors;

%NOTE: The state vector has the following form: x_n = [positionX,
%positionY, velocityX, velocityY]

%TODO: Implement the following hard-coded stuff to be read from parset!!!!!

parset.d0 = parset.simAreaSize/200; %reference distance
% parset.d0 = 1; %reference distance
%NOTE: If the target is closer to a sensor than d0, the measured amplitude 
%is clipped to A/d0. (A is the amplitude in the distance 1 from target)
%This also means that, at distance d0 from the target, the measured 
%amplitude of its sound is A/d0 (and also for smaller distances => this
%prevents the measurements from groving to infinity).

%State transition matrices: (x_n = Phi*x_{n-1} + Gamma*w_n)
parset.Phi = [1 0 1 0;0 1 0 1;0 0 1 0;0 0 0 1];
parset.Gamma = [0.5 0;0 0.5;1 0;0 1];
parset.Gamma0 = [1 0;0 1;0 0;0 0];

%Initial target position and velocity: (also the mean of the prior f(x0))
%NOTE: Currently only up to 4 targets are supported in this initialization
%code.
x0_tmp = [4 4 0.05 0.05;36 36 -0.05 -0.05;165 35 -0.6 0.1;5 195 4 -4].';
if parset.numTargets > 4
    parset.numTargets = 4;
    disp('Unsupported number of targets! Setting the number of targets to 4...')
end
parset.x0 = x0_tmp(:,1:parset.numTargets);
clear x0_tmp;

% %Square root of the covariance of the Gaussian prior f(x0):
% parset.sqrtC0 = sqrtm([0.001 0 0 0;0 0.001 0 0;0 0 0.0000001 0; 0 0 0 0.0000001]);
% % parset.sqrtC0(3:4,:) = zeros(2,4);
parset.sqrtC_w0 = sqrtm(2*[1 0;0 1]);

%Load target trajectory:
if parset.pregenTrajectory == 1 %pre-generated target trajectory
    load x_n_pregen;
    parset.x_n = x_n_pregen;
    clear x_n_pregen;
end

%Load target trajectories for each avg cycle:
if parset.pregenTrajectory == 2 || parset.pregenTrajectory == 3 %randomly pre-generated trajectories for each averaging cycle
    load x_n_rndFixed;
    parset.x_n_rndFixed = x_n_rndFixed;
    clear x_n_rndFixed;
end

if parset.pregenTrajectory == 3 %randomly pre-generated measurements
    load measurements;
    parset.z = measurements;
    clear measurements;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Start simulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simulate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
