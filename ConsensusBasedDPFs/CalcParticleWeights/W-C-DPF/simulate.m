function simulate
%Main simulation function. Responsible for: state generation, measurement
%generation, calling consensus and particle filtering functions, MSE
%calculation and results collection and storing.

global parset;

trajectorySuccessfullyGenerated = 1;

sqrtCu = sqrtm(parset.Cu); %square root of driving noise covariance (used for random number generation)

mse = zeros(parset.avgCycles,parset.simTime,parset.numTargets); %MSE (averaged over all sensors) for each avg. cycle (rows) and for each time instant (cols) and for each target (3rd dimension)
mse_v = zeros(parset.avgCycles,parset.simTime,parset.numTargets); %MSE (velocity) (averaged over all sensors) for each avg. cycle (rows) and for each time instant (cols) and for each target (3rd dimension)

if parset.visualize == 1
    x_n_hat_viz = zeros(parset.numSensors,4*parset.numTargets,parset.simTime);
end

m = 1; %index that counts the number of averaging cycles

while m <= parset.avgCycles
    
    try
      
    clear particleFiltering; %clear persistent variables in function particleFiltering
    
    if parset.dynGenerateTopology == 1
        %generate new network topology for each averaging cycle:
        generateTopology_dyn;%dynamical generation of network topology    
    end
        
    %Due to try-catch it's neccessary to have temporary storage:
    mse_tmp = zeros(parset.numSensors,parset.simTime,parset.numTargets); %MSE (instantaneous error) for each sensor (rows) for each time instant (cols) and for each target (3rd dimension)
    mse_tmp_v = zeros(parset.numSensors,parset.simTime,parset.numTargets); %MSE (velocity) (instantaneous error) for each sensor (rows) for each time instant (cols) and for each target (3rd dimension)
    
    if parset.pregenTrajectory == 0 %randomly generated target trajectory
                
        parset.x_n = zeros(parset.numTargets*4,parset.simTime+1);
        %Generate initial target position:
        for tt = 0:parset.numTargets-1
%             parset.x_n(tt*4+1:(tt+1)*4,1) = parset.x0(:,tt+1) + parset.sqrtC0*randn([4,1]); %the prior f(x0) is a Gaussian with mean parset.x0 and sqrt of covariance equal to parset.sqrtC0
            parset.x_n(tt*4+1:(tt+1)*4,1) = parset.Phi*parset.x0(:,tt+1) + parset.Gamma0*parset.sqrtC_w0*randn([2,1]);
            %Generate new states for each time instant:
            for nn = 2:(parset.simTime+1)
                parset.x_n(tt*4+1:(tt+1)*4,nn) = parset.Phi*parset.x_n(tt*4+1:(tt+1)*4,nn-1) + parset.Gamma*(sqrtCu*randn([2,1]));
                if parset.x_n(tt*4+1,nn) < 0.05*parset.simAreaSize || parset.x_n(tt*4+1,nn) > 0.95*parset.simAreaSize || parset.x_n(tt*4+2,nn) < 0.05*parset.simAreaSize || parset.x_n(tt*4+2,nn) > 0.95*parset.simAreaSize
                    trajectorySuccessfullyGenerated = 0;
                    error('Target is out of the simulation area!!!!');
                end
            end
        end
        
        parset.x_n = parset.x_n(:,2:end); %we do not need the initial state in what follows...
        
        %Enforce minimum distance between targets:------------------------>
        %NOTE: Currently implemented only for 2 targets!!!!
        if parset.numTargets == 2
            for nn = 1:parset.simTime
                if norm(parset.x_n(1:2,nn) - parset.x_n(5:6,nn)) < parset.simAreaSize/4
                    trajectorySuccessfullyGenerated = 0;
                    error('Targets are too close to each other!!!!');                    
                end
            end
        end
        %<-----------------------------------------------------------------        
    end
    
    %Else, if pre-generated target trajectory is used, we do not need to 
    %generate the trajectory, we just use the one stored in parset from the
    %beginning.
    if parset.pregenTrajectory == 2 || parset.pregenTrajectory == 3 %randomly pre-generated trajectories for each averaging cycle
        parset.x_n = parset.x_n_rndFixed(:,:,m); %we set the trajectory for the current averaging cycle
    end
    
    disp(m);
    
    for n = 1:parset.simTime %time recursion
        
        %State transition:------------------------------------------------>
        x_n = parset.x_n(:,n);
        %<-----------------------------------------------------------------       
        
        %Generate measurements:------------------------------------------->
        if parset.pregenTrajectory == 3 %randomly pre-generated measurements
            z_n = parset.z(:,n,m);
        else
            z_n = generateMeasurements(x_n(repmat([1 1 0 0],1,parset.numTargets)==1)); %measurements for each sensor (column vector)
        end
        %<-----------------------------------------------------------------
        
        %Particle filtering for each sensor:------------------------------>
        x_n_hat = particleFiltering(z_n); %estimated target position for each sensor (stored in rows)
        %<-----------------------------------------------------------------
    
        %Results processing:---------------------------------------------->
        for tt = 0:parset.numTargets-1
            tmp = x_n_hat(:,tt*4+1:tt*4+2) - repmat(x_n(tt*4+1:tt*4+2).',parset.numSensors,1);
            mse_tmp(:,n,tt+1) = hypot(tmp(:,1),tmp(:,2)).^2;
            %velocity:
            tmp = x_n_hat(:,tt*4+3:tt*4+4) - repmat(x_n(tt*4+3:tt*4+4).',parset.numSensors,1);
            mse_tmp_v(:,n,tt+1) = hypot(tmp(:,1),tmp(:,2)).^2;            
        end
        %<-----------------------------------------------------------------

        if parset.visualize == 1
            x_n_hat_viz(:,:,n) = x_n_hat;        
        end        
        
    end
    
    mse(m,:,:) = mean(mse_tmp,1);
    mse_v(m,:,:) = mean(mse_tmp_v,1);
    
    m = m + 1;
   
    catch ME
        if trajectorySuccessfullyGenerated == 1 %This ensures that we do not display messages regarding target trajectory generation....
            disp('WARNING in simulate: Exception thrown. Repeating the current averaging cycle.');
            disp(ME.message);
            for st = 1:length(ME.stack)
               funcLineMessage = ['In ' ME.stack(st).name ' on line ' num2str(ME.stack(st).line)];
               disp(funcLineMessage);
            end
            %NOTE: Since m=m+1 is inside the try block, in case of an exception
            %we repeat the current averaging cycle since m wasn't incremented.
        else %trajectorySuccessfullyGenerated == 0
            trajectorySuccessfullyGenerated = 1;
        end
    end
    
    %Storage of partial results:
    if mod(m,floor(parset.avgCycles*0.1)) == 0 %every 10% of avg. cycles store current results
%        mseT = mse / (m - 1);
       mseT = mse;
       save 'mseT' mseT;
       mseT_v = mse_v;
       save 'mseT_v' mseT_v;
    end     
end

save 'mse' mse;
save 'mse_v' mse_v;


if (isunix)
    system('rm mseT.mat'); %delete temporary results to save disk space and speed up copying over network
    system('rm mseT_v.mat');
else
    system('del mseT.mat');
    system('del mseT_v.mat');
end

if parset.visualize == 1
    save 'x_n_hat_viz' x_n_hat_viz;
    x_n_viz = parset.x_n;
    save 'x_n_viz' x_n_viz;
    
    sensorsXY = parset.sensorsPos;
    adjMtx = full(parset.adjacencyMat);
        
    save 'sensorsXY' sensorsXY;
    save 'adjMtx' adjMtx;    
end

end