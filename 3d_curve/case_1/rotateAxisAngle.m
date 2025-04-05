function u = rotateAxisAngle(v, t, angle)
%   This function computes the parallel transport of a tangent vector u from 
%   direction t1 to direction t2 on the unit sphere. If t1 and t2 are aligned,
%   the vector remains unchanged. Otherwise, it is rotated accordingly to 
%   preserve its orientation relative to the surface.
%
%   INPUTS:
%     u  - 3x1 vector to be transported
%     t1 - 3x1 unit tangent vector representing the starting direction
%     t2 - 3x1 unit tangent vector representing the target direction
%
%   OUTPUT:
%     d  - 3x1 vector representing the parallel transport of u from t1 to t2


if (angle == 0)
    u = v;
else
    cs = cos(angle);
	ss = sin(angle);
	u = cs*v + ss*cross(t,v)+dot(t,v)*(1.0-cs)*t;
end
    
end
