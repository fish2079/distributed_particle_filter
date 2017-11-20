try %this deals with the case that neither mse or mseT is found

clear all;

numThresholds = 100;

%load results :----------------------------------------------------------->
try
    load mse;
catch  %#ok<CTCH>
    load mseT;
    numAvgCycles = sum(any(mseT(:,:,1),2))      %#ok<NOPTS>
    mse = mseT(any(mseT(:,:,1),2),:,:); %remove zero rows
end
%<-------------------------------------------------------------------------

[r c w] = size(mse);
rmse = sqrt(reshape(mse,r*c*w,1));

thresholds = linspace(min(rmse),max(rmse),numThresholds);
ccdf = zeros(numThresholds,1);

for t = 1:numThresholds
    ccdf(t) = sum(rmse >= thresholds(t));
end

ccdf = ccdf/(r*c*w);

semilogy(thresholds,ccdf,'r');

save 'ccdf' ccdf;
save 'thresholds' thresholds;

catch %#ok<CTCH>
end