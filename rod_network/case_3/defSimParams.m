function simParams = defSimParams()

simParams = struct();

simParams.totalTime = 5.0;

simParams.dt = 1e-3;

simParams.tol = 1e-3;

simParams.maximum_iter = 10;

simParams.Nsteps = round(simParams.totalTime/simParams.dt);

simParams.plotStep = 100;

end