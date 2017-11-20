function visualizeRMSE_aux(lineColor)


%load results :----------------------------------------------------------->
try
    mse_loc = load('mse.mat');
    mse_loc = mse_loc.mse;
catch  %#ok<CTCH>
    mse_loc = load('mseT.mat');
    mse_loc = mse_loc.mseT;
    numAvgCycles = sum(any(mse_loc(:,:,1),2)) %#ok<NOPRT,NASGU>
    mse_loc = mse_loc(any(mse_loc(:,:,1),2),:,:); %remove zero rows
end
%<-------------------------------------------------------------------------

time = 1:size(mse_loc,2);

targets = size(mse_loc,3);

rmse = sqrt(mean(mse_loc,1));

switch lineColor
    case 1
        targetColors = {'r';'r:';'r--';'r.-'};
    case 2
        targetColors = {'g';'g:';'g--';'g.-'};
    case 3
        targetColors = {'b';'b:';'b--';'b.-'};        
    case 4
        targetColors = {'c';'c:';'c--';'c.-'};
    case 5
        targetColors = {'m';'m:';'m--';'m.-'};
    case 6
        targetColors = {'y';'y:';'y--';'y.-'};
    case 7
        targetColors = {'k';'k:';'k--';'k.-'};
end

figure(21);
hold on;
for tt = 1:targets
    plot(time,rmse(1,:,tt),targetColors{tt});
end

end