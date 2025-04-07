function consParams = defConsParams(plateParams)
% This function defines the constrained parameters based on the physical
% struct's boundary conditions. 
%   Input:
%       plateParams - define a plate struct
%
%   Output:
%       consParams - the constrained parameters contains the constrained
%       index and unconstrained index of the discrete plate model

% Define fixed DOF
fixNode = [17;644];
temp = 1;

for i = 1:2
    for j = 1:3
         consInd(temp) = 3 * (fixNode(i) - 1) + j;
         temp = temp + 1;
    end
   
end

dummyInd = 1:plateParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end
