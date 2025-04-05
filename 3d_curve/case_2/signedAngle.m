function angle = signedAngle( u, v, n )
% This function computes the signed angle between two vectors with respect to a reference normal.
%
%   angle = SIGNEDANGLE(u, v, n) returns the signed angle (in radians) from
%   vector u to vector v, measured in the plane defined by the normal vector n.
%
%   Inputs:
%       u - 3x1 vector (starting direction)
%       v - 3x1 vector (ending direction)
%       n - 3x1 unit normal vector defining the orientation of the plane
%
%   Output:
%       angle - Signed angle in radians, positive if the rotation from u to v
%               follows the right-hand rule with respect to n, negative otherwise.

w = cross(u,v);
angle = atan2( norm(w), dot(u,v) );

if (dot(n,w) < 0) 
    angle = -angle;
end

end