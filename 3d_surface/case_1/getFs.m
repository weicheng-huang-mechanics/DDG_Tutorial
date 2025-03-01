function [Fs, Js] = getFs(plateParams, sElement)

Fs = zeros(plateParams.ndof, 1);
Js = zeros(plateParams.ndof, plateParams.ndof);

for c=1:plateParams.ne
    
    node_1 = getVertex(plateParams.x, sElement(c).nodeIndex(1));
    node_2 = getVertex(plateParams.x, sElement(c).nodeIndex(2));
    
    [dF, dJ] = stretchingForce(node_1, node_2, sElement(c).refLen, sElement(c).EA_local);
    
    index = sElement(c).globalIndex;
    
    Fs(index) = Fs(index) - dF;
    Js(index,index) = Js(index,index) - dJ;
end

end
