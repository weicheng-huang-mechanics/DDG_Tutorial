function [Fg, Jg] = getFg(rodParams)

Fg = rodParams.m .* rodParams.garr;
Jg = zeros(rodParams.ndof, rodParams.ndof);

end
