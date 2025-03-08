function [Fg, Jg] = getFg(rodParams)

garr = zeros(rodParams.ndof, 1);
for c = 1:rodParams.nv
    garr( 2 * (c-1) + 1 : 2 * (c-1) + 2) = rodParams.g;
end

Fg = garr;
Jg = zeros(rodParams.ndof, rodParams.ndof);

end
