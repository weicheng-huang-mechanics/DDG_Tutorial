function angle = signedAngle( u, v, n )
%   This function returns the signed angle between u and v using the right-hand rule 
%   about the normal vector n. The sign of the angle is determined by the direction 
%   of the cross product of u and v relative to n.
%
%   Input:
%     u - 3x1 vector (first direction)
%     v - 3x1 vector (second direction)
%     n - 3x1 normal vector defining the orientation (for determining sign)
%
%   Output:
%     angle - Signed angle (in radians) from u to v, measured around normal n
%

w = cross(u,v);
angle = atan2( norm(w), dot(u,v) );

if (dot(n,w) < 0) 
    angle = -angle;
end

end
