function plotRMSE(filename)
%   Function to plot the boxplot of RMSE of position estimate and run time
%
%   Input:
%       filename: a cell array of strings, each cell contains one filename
%       Each file contains the tracking results for one specific simulation
%       setting (i.e., different values of N)

% initialize output variable
averageError = [];
averageTime = [];

% group variables used for boxplot functions
colorgroup = [];

% xlabel parameters for the plot
xticklabel = {};
xtick = [];
groupSeparator = [];

% Loop through all results
for i=1:numel(filename)
    load(filename{i});
    % Compute RMSE over time
    % The RMSE is stored as Nb_alg x Nb_trials matrix 
    % We concatenante the results from all files
    averageError = [averageError; squeeze(mean(results.pos_error,2))];
    
    % Compute average runtime
    % The runtime is stored as Nb_alg x Nb_trials matrix 
    % We concatenante the results from all files 
    averageTime = [averageTime; results.runtime];
    
    % Set up the grouping variables used in the boxplot functions
    colorgroup = [colorgroup, 1:size(results.runtime,1)];
    
    % Modify as needed depending on actual values
    xticklabel{i} = num2str(parameters.F.N);
    xtick(i) = (size(results.runtime,1)/2+0.5) + (i-1)*size(results.runtime,1);
    groupSeparator(i) = size(results.runtime,1)*i+0.5;
end

% Boxplot of time-averaged RMSE
figure();
set(gcf,'color','white');
set(gca,'fontsize',32);
boxplot(averageError','widths',0.75,'colorgroup',colorgroup,'positions',1:numel(colorgroup));
title('Red - BS, Brown - CSS, Green - LC, Blue - LA, Purple - Cluster');

hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,15],'k');
end

ylim([0,10]);
xlabel('N');
ylabel('RMSE (km)');
% ylabel('Time (s)');

set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);


% Boxplot of total runtime
figure();
set(gcf,'color','white');
set(gca,'fontsize',32);
boxplot(averageTime','widths',0.75,'colorgroup',colorgroup,'positions',1:numel(colorgroup));
title('Red - BS, Brown - CSS, Green - LC, Blue - LA, Purple - Cluster');

hold on;
for i=1:numel(groupSeparator)-1
    plot([groupSeparator(i),groupSeparator(i)],[0,20],'k');
end

ylim([0,20]);
xlabel('N');
ylabel('Time (s)');

set(gca,'xtick', xtick);
set(gca,'xticklabel', xticklabel);