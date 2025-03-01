function simParams = defSimParams()

simParams = struct();

simParams.totalTime = 10.0;

simParams.dt = 1e-2;

simParams.tol = 1e-3;

simParams.maximum_iter = 10;

simParams.Nsteps = round(simParams.totalTime/simParams.dt);

simParams.plotStep = 10;

end