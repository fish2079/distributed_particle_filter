function plotEstimate(S, x_estimate, plotAdjacency)
%   Function to plot the true target track, estimated track and sensor
%   positions
%
%   Input:
%       S: struct containing sensor positions and true track
%       x_estimate: estimated target track, each column represents one time
%       step, first two rows represent estimated position
%       plotAdjacency: boolean flag, wether to plot the sensor neighborhood
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
    h1=plot(S.x_t(1,:), S.x_t(2,:),'b','linewidth',linewidth);
    hold on;
    h2=scatter(S.sensorPos(1,:), S.sensorPos(2,:), scatterSize,'kd','filled');
catch ME
end

if (plotAdjacency)
    for i=1:S.nb_sensors
        for j=i+1:S.nb_sensors
            if (S.A(i,j)==1)
                plot(S.sensorPos(1,[i,j]), S.sensorPos(2,[i,j]),'r--','linewidth',linewidth/2);
            end
        end
    end
end

try
    h3=plot(x_estimate(1,:), x_estimate(2,:), 'b--','linewidth',linewidth);
    legend([h1,h2,h3],{'True track', 'Sensor', 'Estimated track'});
catch ME
    legend([h1,h2],{'True track', 'Sensor'});
end

xlabel('x');
ylabel('y');
axis equal
set(gca,'fontsize',32);