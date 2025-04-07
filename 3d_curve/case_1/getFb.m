function [Fb, Jb] = getFb(rodParams, bElement)
% This function computes the bending force and jacobian of the simulated 
% system. In this case the simulated system is a 3D helix
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%         bElement - the bending element list
%
% OUTPUTS: Fb - bending forces (ndof x 1)
%          Jb - bending jacobian (ndof x ndof)

Fb = zeros(rodParams.ndof, 1);
Jb = zeros(rodParams.ndof, rodParams.ndof);

for c=1:rodParams.nb
    [dF, dJ] = bendingForce(bElement(c).nodePos_1', bElement(c).nodePos_2', bElement(c).nodePos_3', ...
        bElement(c).m_11', bElement(c).m_12', bElement(c).m_21',bElement(c).m_22', bElement(c).kappaBar, ... 
        bElement(c).voroLen, bElement(c).EI_local, bElement(c).directSign_1, bElement(c).directSign_2);
    
    index = bElement(c).globalIndex;
    
    Fb(index) = Fb(index) - dF;
    Jb(index,index) = Jb(index,index) - dJ;
end


end
