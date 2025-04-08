function simParams = defSimParams()

% This function defines the struct that contains the numerical parameters
%   Input:
%
%   Output:
%       simParams - the defined struct contains the numerical parameters 
%                   of the simulated system

simParams = struct();

% Total time
simParams.totalTime = 2.0;

% Time step size
simParams.dt = 1e-2;

% Numerical tolerance
simParams.tol = 1e-4;

% Max iterations 
simParams.maximum_iter = 10;

% Total step
simParams.Nsteps = round(simParams.totalTime/simParams.dt);

% Plot frequency
simParams.plotStep = 1;

end
