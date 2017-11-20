
trackLossTrshld = 5; %[m]

%load results :----------------------------------------------------------->
try
    load mse;
catch  %#ok<CTCH>
    load mseT;
    numAvgCycles = sum(any(mseT(:,:,1),2))      %#ok<NOPTS>    
    mse = mseT(any(mseT(:,:,1),2),:,:); %remove zero rows
end
%<-------------------------------------------------------------------------

save 'mse_ORIG' mse;

[r c w] = size(mse);

mse_max = max(mse,[],3); %choose the max over all targets
indices = sqrt(max(mse_max,[],2))<trackLossTrshld;
% mse_max = max(mse(:,end,:),[],3); %choose the max (at final time) over all targets
% indices = sqrt(mse_max)<trackLossTrshld;

trackLostPercentage = sum(1-indices)/r %#ok<NOPTS>

mse = mse(indices,:,:);

save 'mse' mse;
save 'trackLostPercentage' trackLostPercentage;