function consParams = defConsParams(rodParams)
% This function defines the constrained parameters based on the physical
% struct's boundary conditions. 
%   Input:
%       rodParams - define a cable struct
%
%   Output:
%       consParams - the constrained parameters contains the constrained
%       index and unconstrained index of the discrete beam model

% Define fixed DOF
consInd = [1;2;79;80];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end
