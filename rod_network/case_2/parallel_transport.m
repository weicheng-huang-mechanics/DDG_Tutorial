function d = parallel_transport(u, t1, t2)
%   This function computes the parallel transport of vector u along a 
%   surface or curve, from the local frame defined by tangent vector t1 
%   to the one defined by t2. It preserves the directional relationship 
%   of u relative to t1 as it's mapped to the frame of t2.
%   Inputs:
%       u  - 3x1 vector to be transported
%       t1 - 3x1 unit tangent vector at the starting point
%       t2 - 3x1 unit tangent vector at the ending point
%
%   Output:
%       d  - 3x1 vector after parallel transport of u from t1 to t2

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
