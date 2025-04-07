function d = parallel_transport(u, t1, t2)

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

b = cross(t1, t2);
if (norm(b) == 0 ) 
    d = u;
else
    b = b / norm(b);
    
    b = b - dot(b,t1) * t1;
    b = b / norm(b);
    b = b - dot(b,t2) * t2;
    b = b / norm(b);
    
    n1 = cross(t1, b);
    n2 = cross(t2, b);
    d = dot(u,t1) * t2 + dot(u, n1) * n2 + dot(u, b) * b;
end
end
