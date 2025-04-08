function consParams = defConsParams(rodParams)
% This function defines the constrained parameters based on the physical
% struct's boundary conditions. 
%   Input:
%       rodParams - the defined rod struct contains the physical and
%       numerical parameters of the simulated system
%
%   Output:
%       consParams - the constrained parameters contain the constrained
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
