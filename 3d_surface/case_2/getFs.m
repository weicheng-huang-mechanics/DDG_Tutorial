function [Fs, Js] = getFs(plateParams, sElement)
% This function computes the stretching force and jacobian of the simulated system.
% INPUTS: rodParams - the defined plate struct contains the physical and
%                     numerical parameters of the simulated system
%         sElement - the stretching element list
%
% OUTPUTS: Fs - stretching forces (ndof x 1)
%          Js - stretching jacobian (ndof x ndof)

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
