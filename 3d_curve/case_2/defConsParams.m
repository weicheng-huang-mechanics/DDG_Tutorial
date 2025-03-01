function consParams = defConsParams(rodParams)

nv = rodParams.nv;
ne = rodParams.ne;

% Define fixed DOF
consInd = [1;2;3;4;5;6;3*nv;3*nv-1;3*nv-2;3*nv-3;3*nv-4;3*nv+1;3*nv+ne];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end