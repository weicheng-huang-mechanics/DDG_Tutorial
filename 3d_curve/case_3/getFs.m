function [Fs, Js] = getFs(rodParams, sElement)

Fs = zeros(rodParams.ndof, 1);
Js = zeros(rodParams.ndof, rodParams.ndof);

for c=1:rodParams.ne
    
    [dF, dJ] = stretchingForce(sElement(c).nodePos_1, sElement(c).nodePos_2, sElement(c).refLen, sElement(c).EA_local);
    
    index = sElement(c).globalIndex;
    
    Fs(index) = Fs(index) - dF;
    Js(index,index) = Js(index,index) - dJ;
end

end
