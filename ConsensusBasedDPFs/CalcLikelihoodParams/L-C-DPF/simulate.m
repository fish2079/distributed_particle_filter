function simulate
%Main simulation function. Responsible for: state generation, measurement
%generation, calling consensus and particle filtering functions, MSE
%calculation and results collection and storing.

global parset;
numberOfLostTracks = 0; %used to count lost tracks if track loss detection is used

trajectorySuccessfullyGenerated = 1;

sqrtCu = sqrtm(parset.Cu); %square root of driving noise covariance (used for random number generation)

% mse = zeros(parset.numSensors,parset.simTime); %MSE for each sensor (rows) for each time instant (cols)

% mse = zeros(parset.numSensors,parset.simTime,parset.avgCycles); %MSE for each sensor (rows) for each time instant (cols) and for each avg. cycle (3rd dimension)
% %NOTE: The entries of the matrix above are not really MSE (since they are
% %not averaged); they are just instantaneous errors.

mse = zeros(parset.avgCycles,parset.simTime,parset.numTargets); %MSE (averaged over all sensors) for each avg. cycle (rows) and for each time instant (cols) and for each target (3rd dimension)
mse_v = zeros(parset.avgCycles,parset.simTime,parset.numTargets); %MSE (velocity) (averaged over all sensors) for each avg. cycle (rows) and for each time instant (cols) and for each target (3rd dimension)

if parset.visualize == 1
    x_n_hat_viz = zeros(parset.numSensors,4*parset.numTargets,parset.simTime);
end

m = 1; %index that counts the number of averaging cycles

while m <= parset.avgCycles
    
    try
        
    switch parset.pfAlg
        case 0 %SIR
            clear particleFiltering_SIR; %clear persistent variables
        case 1 %GPF
            clear particleFiltering_GPF; %clear persistent variables
        case 2 %GSPF
            clear particleFiltering_GSPF; %clear persistent variables
        otherwise
            error('Wrong setting of Particle filter algorithm');
    end    

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
        
        %Particle filter prediction step:--------------------------------->
        %NOTE: The predicted particles are stored also in the global 
        %variable particlesApprox and are used for reference point 
        %calculation in Taylor approximation and for least squares 
        %approximation in likelihood consensus.
        switch parset.pfAlg
            case 0 %SIR
                particleFiltering_SIR(1);
            case 1 %GPF
                particleFiltering_GPF(1);
            case 2 %GSPF
                particleFiltering_GSPF(1);
        end
        %<-----------------------------------------------------------------
        
        %Joint likelihood calculation:------------------------------------>
        switch parset.consensusType
            case 0 %no consensus (centralized and exact calculation of joint likelihood)
                jointLikelihoodSums = z_n.'; %contains measurements which are then used to calculate true joint likelihood
            case 1 %perfect likelihood consensus
                jointLikelihoodSums = likelihoodConsensus_perfect(z_n); %sums in joint likelihood for each sensor (stored in rows)
            case 2 %realistic likelihood consensus
                jointLikelihoodSums = likelihoodConsensus_real(z_n); %sums in joint likelihood for each sensor (stored in rows)
        end
        %<-----------------------------------------------------------------
        
        %Particle filter measurement update:------------------------------>
        switch parset.pfAlg
            case 0 %SIR
                x_n_hat = particleFiltering_SIR(0,jointLikelihoodSums); %estimated target position for each sensor (stored in rows)                
            case 1 %GPF
                x_n_hat = particleFiltering_GPF(0,jointLikelihoodSums); %estimated target position for each sensor (stored in rows)                
            case 2 %GSPF
                x_n_hat = particleFiltering_GSPF(0,jointLikelihoodSums); %estimated target position for each sensor (stored in rows)                
        end
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
        
        %Clairvoyant track loss detection:
        if sqrt(max(max(mse_tmp(:,n,:)))) > parset.trackLossDetTrshld && parset.trackLossDetOnOff == 1
            numberOfLostTracks = numberOfLostTracks + 1;
            error('Track lost!');
        end        
        
        if parset.visualize == 1
            x_n_hat_viz(:,:,n) = x_n_hat;        
        end
    end
    
    
%     mse = mse + mse_tmp;

%     mse(:,:,m) = mse_tmp;

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

% mse = mse / (m - 1);

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

if parset.trackLossDetOnOff == 1
    trackLostPercentage = numberOfLostTracks / (parset.avgCycles+numberOfLostTracks);
    save 'trackLostPercentage' trackLostPercentage;    
end

end