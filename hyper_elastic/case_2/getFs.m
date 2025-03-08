function [Fs, Js] = getFs(rodParams, sElement)


Fs = zeros(rodParams.ndof, 1);
Js = zeros(rodParams.ndof, rodParams.ndof);

for c=1:rodParams.ne
    
    node_1 = getVertex(rodParams.x, sElement(c).nodeIndex(1));
    node_2 = getVertex(rodParams.x, sElement(c).nodeIndex(2));
    
    [dF, dJ] = stretchingForce(node_1, node_2, sElement(c).refLen, sElement(c).refRadius, sElement(c).refThickness, 0.4 * sElement(c).E_local, 0.1 * sElement(c).E_local);
    
    
    index = sElement(c).globalIndex;
    
    Fs(index) = Fs(index) - dF;
    Js(index,index) = Js(index,index) - dJ;
end

end
