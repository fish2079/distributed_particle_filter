% This file is necessary for running the simulations remotely. On the
% remote host type: 
% nohup matlab -nodisplay < runRemotely.m > output.txt &
%
% Important: when starting particular m-file directly from shell, it is
% important that the m-file is a script and not a function.

main; % just call the main...

exit; %exit MATLAB -> it is particularily usefull on Mac OS servers since on them MATLAB usualy doesn't terminate after the job is finished