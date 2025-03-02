function [Fb, Jb] = getFb(rodParams, bElement)

Fb = zeros(rodParams.ndof, 1);
Jb = zeros(rodParams.ndof, rodParams.ndof);

for c=1:rodParams.nb
    
    node_1 = zeros(3,1);
    node_2 = zeros(3,1);
    node_3 = zeros(3,1);

    node_1(1:2) = getVertex(rodParams.x, bElement(c).nodeIndex(1));
    node_2(1:2) = getVertex(rodParams.x, bElement(c).nodeIndex(2));
    node_3(1:2) = getVertex(rodParams.x, bElement(c).nodeIndex(3));
    
    node_1(3) = 0.0;
    node_2(3) = 0.0;
    node_3(3) = 0.0;
    
    [dF, dJ] = bendingForce(node_1, node_2, node_3, bElement(c).kappaBar_1, bElement(c).kappaBar_2, bElement(c).refRadius, bElement(c).voroLen, bElement(c).EI_local);
    
    index = bElement(c).globalIndex;
    
    Fb(index) = Fb(index) - dF;
    Jb(index,index) = Jb(index,index) - dJ;
end


end
