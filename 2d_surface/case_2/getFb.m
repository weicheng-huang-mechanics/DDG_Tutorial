function [Fb, Jb] = getFb(rodParams, bElement)

% This function computes the bending force and jacobian of the simulated system.
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%         bElement - the bending element list
%
% OUTPUTS: Fb - bending forces (2*nv x 1)
%          Jb - bending jacobian (2*nv x 2*nv)

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
