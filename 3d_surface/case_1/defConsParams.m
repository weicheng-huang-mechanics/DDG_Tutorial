function consParams = defConsParams(plateParams)

consParams = struct();

% Define fixed DOF
consInd1 = zeros(3,1);
temp = 1;
temp2 = 1;
for i = 1:plateParams.nv
    xCurrent = getVertex(plateParams.x, i);
    
    if (xCurrent(1) < 0.15)
        consInd1(temp) = i;
        temp = temp + 1;
        
        consInd(temp2) = 3 * (i-1) + 1;
        temp2 = temp2 + 1;
        consInd(temp2) = 3 * (i-1) + 2;
        temp2 = temp2 + 1;
        consInd(temp2) = 3 * (i-1) + 3;
        temp2 = temp2 + 1;
    end
end

% consInd2 = zeros(3,1);
% temp = 1;
% for i = 1:plateParams.nv
%     xCurrent = getVertex(plateParams.x, i);
%     
%     if (xCurrent(1) > 0.94)
%         consInd2(temp) = i;
%         temp = temp + 1;
%         
%         consInd(temp2) = 3 * (i-1) + 1;
%         temp2 = temp2 + 1;
%         consInd(temp2) = 3 * (i-1) + 2;
%         temp2 = temp2 + 1;
%         consInd(temp2) = 3 * (i-1) + 3;
%         temp2 = temp2 + 1;
%     end
% end

consParams.consInd1 = consInd1;
%consParams.consInd2 = consInd2;

[consParams.fixNv,~] = size(consParams.consInd1);


dummyInd = 1:plateParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end