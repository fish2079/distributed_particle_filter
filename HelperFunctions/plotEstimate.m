function plotEstimate(S, x_estimate)
%   Function to plot the true target track, estimated track and sensor
%   positions
%
%   Input:
%       S: struct containing sensor positions and true track
%       x_estimate: estimated target track, each column represents one time
%       step, first two rows represent estimated position
%
% Jun Ye Yu
% McGill University
% jun.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

linewidth = 6;
scatterSize = 500;
figure();
set(gcf,'color','white');
try
    plot(S.x_t(1,:), S.x_t(2,:),'b','linewidth',linewidth);
    hold on;
    scatter(S.sensorPos(1,:), S.sensorPos(2,:), scatterSize,'kd','filled');
catch ME
end
try
    plot(x_estimate(1,:), x_estimate(2,:), 'r--','linewidth',linewidth);
    legend('True track', 'Sensor', 'Estimated track');
catch ME
    legend('True track', 'Sensor');
end



xlabel('x');
ylabel('y');
set(gca,'fontsize',32);