function consParams = defConsParams(plateParams)

% Define fixed DOF
%fixNode = [10;11;21;252;262;263];
fixNode = [17;644];
temp = 1;

for i = 1:2
    for j = 1:3
         consInd(temp) = 3 * (fixNode(i) - 1) + j;
         temp = temp + 1;
    end
   
end
%consInd = [10*3+1;10*3+2;10*3+3;262*3+1;262*3+2;262*3+3];

dummyInd = 1:plateParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end