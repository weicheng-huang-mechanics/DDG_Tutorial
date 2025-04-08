function consParams = defConsParams(rodParams)

% This function defines the constrained parameters based on the physical
% struct's boundary conditions. 
%   Input:
%       rodParams - define a shell struct
%
%   Output:
%       consParams - the constrained parameters contain the constrained
%       index and unconstrained index of the discrete shell mode

% Define fixed DOF
consInd = [1;2;3;4;2*rodParams.nv];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end
