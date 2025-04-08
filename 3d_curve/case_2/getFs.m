function [Fs, Js] = getFs(rodParams, sElement)
% This function computes the stretching force and jacobian of the simulated system. 
% Input:  rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%         sElement - the stretching element list
%
% Output:  Fs - stretching forces (ndof x 1)
%          Js - stretching jacobian (ndof x ndof)

Fs = zeros(rodParams.ndof, 1);
Js = zeros(rodParams.ndof, rodParams.ndof);

for c=1:rodParams.ne
    
    [dF, dJ] = stretchingForce(sElement(c).nodePos_1, sElement(c).nodePos_2, sElement(c).refLen, sElement(c).EA_local);
    
    index = sElement(c).globalIndex;
    
    Fs(index) = Fs(index) - dF;
    Js(index,index) = Js(index,index) - dJ;
end

end
