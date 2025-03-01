function consParams = defConsParams(rodParams)

%nv = rodParams.nv;
%ne = rodParams.ne;

% Define fixed DOF
consInd = [2;3;31;33;62;63;91;93];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end