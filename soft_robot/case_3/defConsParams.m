function consParams = defConsParams(rodParams)
% This function defines the constrained parameters based on the physical
% struct's boundary conditions. 
%   Input:
%       rodParams - define a beam struct
%
%   Output:
%       consParams - the constrained parameters contain the constrained
%       index and unconstrained index of the discrete beam model

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
