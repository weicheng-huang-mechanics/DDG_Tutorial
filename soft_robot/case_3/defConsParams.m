function consParams = defConsParams(rodParams)

% Define fixed DOF
%consInd = linspace(1,2*49,2*49);

consInd = [];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end