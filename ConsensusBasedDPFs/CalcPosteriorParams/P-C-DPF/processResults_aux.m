clear all;

%load results :----------------------------------------------------------->
try
    load mse;
catch  %#ok<CTCH>
    load mseT;
    numAvgCycles = sum(any(mseT(:,:,1),2))      %#ok<NOPTS>    
    mse = mseT(any(mseT(:,:,1),2),:,:); %remove zero rows
end
%<-------------------------------------------------------------------------

save 'mse_ORIG' mse; %save the original mse file...

[r c w] = size(mse);

for tt = 1:w
    
disp('Target:');
disp(tt);

figure(15);

mse_tt = mse(:,:,tt);

cont = 1;
while cont > 0


%[value, maxAvgCycle] = max(max(mse_tt,[],2)); %find the avg cycle where max MSE is reached
[value, maxAvgCycle] = max(mse_tt(:,end)); %find the avg cycle where max MSE at the final time is reached


subplot(2,1,1);
plot(sqrt(mse_tt(maxAvgCycle,:)));

subplot(2,1,2);
plot(sqrt(mean(mse(:,:,tt))));
hold on
indices = (ones(1,size(mse_tt,1)) == 1);
indices(maxAvgCycle) = (0 == 1);
plot(sqrt(mean(mse_tt(indices,:))),'r');
hold off

cont = input('To stop type negative value. To continue type positive. Input: ');

mse_tt = mse_tt(indices,:);
mse = mse(indices,:,:);

end

end

save 'mse' mse;