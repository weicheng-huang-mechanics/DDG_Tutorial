function [Fs, Js] = getFs(rodParams, sElement)
% This function computes the stretching force and jacobian of the simulated system.
% In this case the simulated system is a rotational membrane
% Input:  rodParams - the defined membrane struct contains the physical and
%                     numerical parameters of the simulated system
%         sElement - the stretching element list (ne x 2)
%
% Output: Fs - stretching forces (ndof x 1)
%         Js - stretching jacobian (ndof x ndof)

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
