function consParams = defConsParams(rodParams)
% This function defines the constrained parameters based on the physical
% struct's boundary conditions. 
%   Input:
%       rodParams - define a rod struct
%
%   Output:
%       consParams - the constrained parameters contains the constrained
%       index and unconstrained index of the discrete beam model

% Define fixed DOF
fixIndex = [1;29;30;72;73;121;122;164;165;193;194;219;220;257;258;301;302;339;340;365];

temp = 1;
consInd = zeros(3,1);
for i = 1:size(fixIndex)
    consInd(temp) = 3 * (fixIndex(i)-1) + 1;
    temp = temp + 1;
    consInd(temp) = 3 * (fixIndex(i)-1) + 2;
    temp = temp + 1;
    consInd(temp) = 3 * (fixIndex(i)-1) + 3;
    temp = temp + 1;
end
% 
% for i = 3*rodParams.nv+1:rodParams.ndof;
%     consInd(temp) = i;
%     temp = temp + 1;
% end

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

consParams.fixIndex = fixIndex;

end