function [Fg, Jg] = getFg(rodParams)
% This function computes the gravitational force and jacobian of the simulated 
% system. In this case the simulated system is a 2D beam
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%
% OUTPUTS: Fg - gravitational forces (2*nv x 1)
%          Jg - gravitational jacobian (2*nv x 2*nv)

garr = zeros(rodParams.ndof, 1);
for c = 1:rodParams.nv
    garr( 2 * (c-1) + 1 : 2 * (c-1) + 2) = rodParams.g;
end

Fg = garr;
Jg = zeros(rodParams.ndof, rodParams.ndof);

end
