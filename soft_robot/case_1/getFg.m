function [Fg, Jg] = getFg(rodParams)
% This function computes the gravitational force and jacobian of the simulated system.
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%
% OUTPUTS: Fg - gravitational forces (2*nv x 1)
%          Jg - gravitational jacobian (2*nv x 2*nv)

Fg = rodParams.m .* rodParams.garr;
Jg = zeros(rodParams.ndof, rodParams.ndof);

end
