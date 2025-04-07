function simParams = defSimParams()
% This function defines the a struct contains the numerical parameters
%   Input:
%
%   Output:
%       simParams - the defined struct contains the numerical parameters 
%                   of the simulated system

simParams = struct();

simParams.totalTime = 10.0;

simParams.dt = 1e-2;

simParams.tol = 1e-3;

simParams.maximum_iter = 10;

simParams.Nsteps = round(simParams.totalTime/simParams.dt);

simParams.plotStep = 10;

end