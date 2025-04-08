function [Fg, Jg] = getFg(rodParams)
% This function computes the gravitational force and jacobian of the simulated system.
% Input:  rodParams - the defined beam struct contains the physical and
%                     numerical parameters of the simulated system
%
% Output: Fg - gravitational forces (ndof x 1)
%         Jg - gravitational jacobian (ndof x ndof)

Fg = rodParams.m .* rodParams.garr;
Jg = zeros(rodParams.ndof, rodParams.ndof);

end
