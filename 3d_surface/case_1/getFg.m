function [Fg, Jg] = getFg(plateParams)
% This function computes the gravitational force and jacobian of the simulated system.
% Input:  plateParams - the defined plate struct contains the physical and
%                     numerical parameters of the simulated system
%
% Output: Fg - gravitational forces (ndof x 1)
%         Jg - gravitational jacobian (ndof x ndof)


Fg = plateParams.m .* plateParams.garr;
Jg = zeros(plateParams.ndof, plateParams.ndof);

end
