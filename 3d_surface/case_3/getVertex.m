function localNode = getVertex(x, i)
% This function extracts the 3d coordinate of the i-th edge from the global 
% DOF vector.
%
%   Input:
%       x - Global degrees of freedom vector (ndof x 1)
%       i - Index of the node to extract (1-based)
%
%   Output:
%       localNode - a scalar containing the i-th node's coordinate

localNode = zeros(3,1);

localNode(1) = x(3*(i-1)+1);
localNode(2) = x(3*(i-1)+2);
localNode(3) = x(3*(i-1)+3);

end
