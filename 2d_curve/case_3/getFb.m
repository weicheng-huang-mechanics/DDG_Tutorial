function [Fb, Jb] = getFb(rodParams, bElement)
% This function computes the bending force and jacobian of the simulated 
% system. In this case the simulated system is a 2D beam
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%         bElement - the bending element list (nb x 2)
%
% OUTPUTS: Fb - bending forces (2*nv x 1)
%          Jb - bending jacobian (2*nv x 2*nv)

Fb = zeros(rodParams.ndof, 1);
Jb = zeros(rodParams.ndof, rodParams.ndof);

mVec = [0;0;1];

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
    
    [dF, dJ] = bendingForce(node_1, node_2, node_3, mVec, mVec, bElement(c).nBar, bElement(c).voroLen, bElement(c).EI_local);
    
    index = bElement(c).globalIndex;
    
    Fb(index) = Fb(index) - dF;
    Jb(index,index) = Jb(index,index) - dJ;
end


end
