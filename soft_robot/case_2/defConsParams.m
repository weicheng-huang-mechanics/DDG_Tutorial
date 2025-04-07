function consParams = defConsParams(rodParams)
% This function defines the constrained parameters based on the physical
% struct's boundary conditions. 
%   Input:
%       rodParams - define a rod struct (beam is a 2D rod)
%
%   Output:
%       consParams - the constrained parameters contains the constrained
%       index and unconstrained index of the discrete beam model

% Define fixed DOF
consInd = [80];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end