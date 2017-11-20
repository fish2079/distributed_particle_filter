function generateTopology_dyn
%This function allows the network topology to change dynamically from
%averaging cycle to averaging cycle. Since it uses global variables, it is
%called at the beginning of each averaging cycle with no arguments.
%
%This function generates the 2D physical topology of the sensor network. It
%changes two matrices in global parset. The first one being sensorsPos 
%(numSensors x 2 matrix) containing XY coordinates of every sensor 
%(in rows), and the second one being adjacencyMat (numSensors x numSensors 
%matrix) defining neighborhood (adjacency matrix of the corresponding 
%graph).
%
%Sensors are placed at "noisy rectangular grid" positions. 

global parset;

%Simulation parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSensors_oneDir = 5; % the number of sensors in one direction (mnumber of sensors in one direction of the noisy rectangular grid)
posDeviation = 2; % controls the amount of "deviation" from the rectangular grid positions
neighborhoodRadius = 18; % the neighborhood of a given sensor is given by those sensors that are nearer than this value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simAreaSize = 40; % the simulation area is a square of simAreaSize x simAreaSize [m]
%NOTE: This value must be the same as the simAreaSize value found in main!

numSensors = numSensors_oneDir * numSensors_oneDir; % total number of sensors

gridStep = simAreaSize / (numSensors_oneDir - 1);

sensorsXY = zeros(numSensors,2);
sensorsXY(:,1) = repmat((0:(numSensors_oneDir-1))*gridStep,1,numSensors_oneDir);
sensorsXY(:,2) = reshape(repmat((0:(numSensors_oneDir-1))*gridStep,numSensors_oneDir,1),1,numSensors);
sensorsXY = sensorsXY + randn(numSensors,2)*posDeviation;

% make sure that sensor x and y coordinates are in the interval <0,simAreaSize> 
sensorsXY(sensorsXY < 0) = 0;
sensorsXY(sensorsXY > simAreaSize) = simAreaSize;

% create adjacency matrix determined by the sensor positions and the neighborhood size:
adjMtx = zeros(numSensors);

for iter_sensor=1:numSensors 
    
    for iter_sensor2=1:numSensors 
        
        if (norm(sensorsXY(iter_sensor,:)-sensorsXY(iter_sensor2,:))<= neighborhoodRadius) && (iter_sensor ~= iter_sensor2)
            
            adjMtx (iter_sensor,iter_sensor2) = 1 ; 
            
        end
    end
end

parset.sensorsPos = sensorsXY; %physical positions of the sensors  
parset.adjacencyMat = sparse(adjMtx);  %adjacency matrix of the graph representation of the sensor network (with communication links)

end
