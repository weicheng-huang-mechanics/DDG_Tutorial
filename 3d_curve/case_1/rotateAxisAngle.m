function u = rotateAxisAngle(v, t, angle)

%   This function computes the rotation of a vector along another vector 
%   by an angle.
%
%   INPUTS:
%     v - 3x1 vector to be rotated
%     t - 3x1 unit tangent vector 
%     angle - a scale for the rotational angle
%
%   OUTPUT:
%     u  - 3x1 vector representing the rotational vector

if (angle == 0)
    u = v;
else
    cs = cos(angle);
	ss = sin(angle);
	u = cs*v + ss*cross(t,v)+dot(t,v)*(1.0-cs)*t;
end
    
end
