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
%fixNode = [];

fixNode = zeros(3,1);

temp = 1;
for i = 1:plateParams.nv

    if (plateParams.x(3*(i-1)+1) < 0.40)
        fixNode(temp) = i;
        temp = temp + 1;
    end

%     if (plateParams.x(3*(i-1)+1) > 2.8)
%         fixNode(temp) = i;
%         temp = temp + 1;
%     end
end

consInd = zeros(3,1);
temp2 = 1;
for i = 1:temp-1
    for j = 1:3
         consInd(temp2) = 3 * (fixNode(i) - 1) + j;
         temp2 = temp2 + 1;
    end
end

consInd = [consInd;(417-1) * 3 + 3];

dummyInd = 1:plateParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end
