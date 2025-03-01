function [Fg, Jg] = getFg(rodParams)

garr = zeros(rodParams.ndof, 1);
for c = 1:rodParams.nv
    garr( 3 * (c-1) + 1 : 3 * (c-1) + 3) = rodParams.g;
end

Fg = rodParams.m .* garr;
Jg = zeros(rodParams.ndof, rodParams.ndof);

end
