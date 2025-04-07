function [Fg, Jg] = getFg(rodParams)
% This function computes the gravitational force and jacobian of the simulated 
% system. In this case the simulated system is a 2D beam
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%
% OUTPUTS: Fg - gravitational forces (ndof x 1)
%          Jg - gravitational jacobian (ndof x ndof)

Fg = rodParams.m .* rodParams.garr;
Jg = zeros(rodParams.ndof, rodParams.ndof);

end
