function consParams = defConsParams(rodParams)

% Define fixed DOF
consInd = [1;2;3;4;2*rodParams.nv-3;2*rodParams.nv-2;2*rodParams.nv-1;2*rodParams.nv];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end