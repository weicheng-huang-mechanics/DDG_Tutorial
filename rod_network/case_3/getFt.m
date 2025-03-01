function [Ft, Jt] = getFt(rodParams, bElement)

Ft = zeros(rodParams.ndof, 1);
Jt = zeros(rodParams.ndof, rodParams.ndof);

for c=1:rodParams.nb
    [dF, dJ] = twistingForce(bElement(c).nodePos_1', bElement(c).nodePos_2', bElement(c).nodePos_3', ...
        bElement(c).theta_1, bElement(c).theta_2, bElement(c).refTwist, bElement(c).voroLen, ... 
        bElement(c).EI_local, bElement(c).directSign_1, bElement(c).directSign_2);
    
    index = bElement(c).globalIndex;
    
    Ft(index) = Ft(index) - dF;
    Jt(index,index) = Jt(index,index) - dJ;
end


end
