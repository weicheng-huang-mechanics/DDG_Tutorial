function localNode = getTheta(x, nv, i)
% This function extracts the rotation angle of the i-th edge from the global 
% DOF vector.
%
%   Input:
%       x - Global degrees of freedom vector (ndof x 1)
%       nv - number of the vertices
%       i - Index of the edge to extract (1-based)
%
%   Output:
%       localNode - a scalar containing the i-th edge's theta


localNode = x(3*nv+i);

end
