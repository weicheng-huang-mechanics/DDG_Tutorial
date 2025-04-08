function localNode = getVertex(x, i)
% This function extracts the 2D position of the i-th node from the global 
% DOF vector.
%
%   Input:
%       x - global degrees of freedom vector (2*nv x 1), where each node 
%           has 2 DOFs
%       i - index of the node to extract (1-based)
%
%   Output:
%       localNode - 2x1 vector containing the i-th node's coordinate

localNode = zeros(2,1);

localNode(1) = x(2*(i-1)+1);
localNode(2) = x(2*(i-1)+2);

end
