function fig_handle = visualizeParticles(S, D, x_estimate, particles, plotBearing, fig_handle)
%   Function to plot the true target track, sensor positions, particle
%   cloud and bearing measurements
%
%   Input:
%       S: struct containing sensor positions and true track
%       D: struct containing sensor measurements
%       x_estimate: estimated track, each column represents one time step
%       particles: particle cloud, each column represents one particle,
%       first two rows represent estimated position, last row is particle
%       weight
%       plotBearing: boolean flag to plot bearing measurements
%
%   Output:
%       fig_handle: hanle of the figure
%
% Jun Ye Yu
% McGill University
% jun.y.yu@mail.mcgill.ca
% Nov. 9th, 2017

k = size(x_estimate,2);
linewidth = 6;
scatterSize = 500;
if(isempty(fig_handle))
    fig_handle = figure();
else
    clf(fig_handle);
end
set(gcf,'color','white');
set(gca,'fontsize',32);

scatter(particles(1,:), particles(2,:), scatterSize/2, particles(5,:)); colorbar;

hold on;

plot(x_estimate(1,1:k), x_estimate(2,1:k), 'r--','linewidth',linewidth);

plot(S.x_t(1,1:k), S.x_t(2,1:k),'b','linewidth',linewidth);

scatter(S.sensorPos(1,:), S.sensorPos(2,:), scatterSize,'kd','filled');

scatter(x_estimate(1,k), x_estimate(2,k), scatterSize, 'rd','filled');

% Plot bearings if necessary
if (plotBearing)
    % Compute range
    range = sqrt(sum((D.sensorLoc-x_estimate(1:2,k)).^2,1));
    quiver(D.sensorLoc(1,:), D.sensorLoc(2,:), sin(D.measurements).*range, cos(D.measurements).*range, 0);
end
legend('Particles', 'Estimated track', 'True track');

xlabel('x');
ylabel('y');