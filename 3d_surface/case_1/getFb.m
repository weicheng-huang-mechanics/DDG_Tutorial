function [Fb, Jb] = getFb(plateParams, bElement)

Fb = zeros(plateParams.ndof, 1);
Jb = zeros(plateParams.ndof, plateParams.ndof);

for c=1:plateParams.nb
    
    node_1 = getVertex(plateParams.x, bElement(c).nodeIndex(1));
    node_2 = getVertex(plateParams.x, bElement(c).nodeIndex(2));
    node_3 = getVertex(plateParams.x, bElement(c).nodeIndex(3));
    node_4 = getVertex(plateParams.x, bElement(c).nodeIndex(4));
    
    [dF, dJ] = bendingForce(node_1, node_2, node_3, node_4, bElement(c).nBar, bElement(c).EI_local);
    
    index = bElement(c).globalIndex;
    
    Fb(index) = Fb(index) + dF;
    Jb(index,index) = Jb(index,index) + dJ;
end


end
