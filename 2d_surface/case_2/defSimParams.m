function simParams = defSimParams()

simParams = struct();

% Total time
simParams.totalTime = 6.0;

% Time step size
simParams.dt = 5e-3;

% Numerical tolerance
simParams.tol = 1e-4;

% Max iterations 
simParams.maximum_iter = 10;

% Total step
simParams.Nsteps = round(simParams.totalTime/simParams.dt);

% Plot frequency
simParams.plotStep = 1;

end
