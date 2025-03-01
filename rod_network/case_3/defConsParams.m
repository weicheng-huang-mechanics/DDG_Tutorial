function consParams = defConsParams(rodParams)

consInd = [1;2;3;4;5;6;31;32;33;58;59;60;3*rodParams.nv+1];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end