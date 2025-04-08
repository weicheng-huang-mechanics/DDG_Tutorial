function [Fb, Jb] = getFb(plateParams, bElement)
% This function computes the bending force and jacobian of the simulated system.
% Input:  plateParams - the defined plate struct contains the physical and
%                     numerical parameters of the simulated system
%         bElement - the bending element list
%
% Output: Fb - bending forces (ndof x 1)
%         Jb - bending jacobian (ndof x ndof)

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
