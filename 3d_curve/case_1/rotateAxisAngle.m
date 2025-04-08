function u = rotateAxisAngle(v, t, angle)
%   This function rotates the input vector v around the 
%   axis t by the given angle (in radians), using the axis-angle rotation formula.
%
%   Input:
%       v     - 3x1 input vector to be rotated
%       t     - 3x1 unit vector representing the axis of rotation
%       angle - scalar rotation angle in radians
%
%   Output:
%       u     - 3x1 rotated vector

if (angle == 0)
    u = v;
else
    cs = cos(angle);
	ss = sin(angle);
	u = cs*v + ss*cross(t,v)+dot(t,v)*(1.0-cs)*t;
end
    
end
