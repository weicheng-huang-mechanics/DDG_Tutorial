function consParams = defConsParams(rodParams)

% Define fixed DOF
fixIndex = [2;21;36;51;66;81];
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

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end