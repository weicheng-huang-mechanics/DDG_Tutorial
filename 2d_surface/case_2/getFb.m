function [Fb, Jb] = getFb(rodParams, bElement)

Fb = zeros(rodParams.ndof, 1);
Jb = zeros(rodParams.ndof, rodParams.ndof);

for c=1:rodParams.nb
    
    node_1 = getVertex(rodParams.x, bElement(c).nodeIndex(1));
    node_2 = getVertex(rodParams.x, bElement(c).nodeIndex(2));
    node_3 = getVertex(rodParams.x, bElement(c).nodeIndex(3));
    
    [dF, dJ] = bendingForce(node_1, node_2, node_3, bElement(c).kappaBar_1, bElement(c).kappaBar_2, bElement(c).voroLen, bElement(c).refRadius, bElement(c).EI_local, bElement(c).nu_local);
    
    index = bElement(c).globalIndex;
    
    Fb(index) = Fb(index) - dF;
    Jb(index,index) = Jb(index,index) - dJ;
end


end
