function d = parallel_transport(u, t1, t2)

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
