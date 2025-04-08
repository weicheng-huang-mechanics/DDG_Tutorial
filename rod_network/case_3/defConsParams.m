function consParams = defConsParams(rodParams)
% This function defines the constrained parameters based on the physical
% struct's boundary conditions. 
%   Input:
%       rodParams - define a rod struct
%
%   Output:
%       consParams - the constrained parameters contain the constrained
%       index and unconstrained index of the discrete rodd model

consInd = [1;2;3;4;5;6;31;32;33;58;59;60;3*rodParams.nv+1];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end
