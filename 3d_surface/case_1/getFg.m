function [Fg, Jg] = getFg(plateParams)

Fg = plateParams.m .* plateParams.garr;
Jg = zeros(plateParams.ndof, plateParams.ndof);

end
