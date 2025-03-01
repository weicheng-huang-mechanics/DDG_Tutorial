function [Fd, Jd] = getFd(rodParams, sElement)


Fd = zeros(rodParams.ndof, 1);
Jd = zeros(rodParams.ndof, rodParams.ndof);

viscosityA = 0.01;
viscosityB = 1.0;


for c=1:rodParams.ne
    
    node_1 = getVertex(rodParams.x, sElement(c).nodeIndex(1));
    node_2 = getVertex(rodParams.x, sElement(c).nodeIndex(2));
    
    edge = norm(node_2 - node_1);
    
    t_current = (node_2 - node_1) / norm(node_2 - node_1);
    n_current = [-t_current(2); t_current(1)];
    
    velo_1 = (getVertex(rodParams.x, sElement(c).nodeIndex(1)) - getVertex(rodParams.x0, sElement(c).nodeIndex(1))) / rodParams.dt;
    velo_2 = (getVertex(rodParams.x, sElement(c).nodeIndex(1)) - getVertex(rodParams.x0, sElement(c).nodeIndex(1))) / rodParams.dt;
    
    index = sElement(c).globalIndex;
    
    velo = velo_1;
    
    uNormal = dot(n_current,velo);
    if (uNormal < 0)
        uNormal   = - uNormal;
        n_current = - n_current;
    end
    
    uTangent = dot(t_current,velo);
    if (uTangent < 0)
        uTangent  = - uTangent;
        t_current = - t_current;
    end
    
    fLocal = - 0.5 * 1000 * edge * (2 * rodParams.r0) * (viscosityA * uTangent * uTangent * t_current + viscosityB * uNormal * uNormal * n_current);
    
    jac = - 1000 * edge * (2 * rodParams.r0) * (viscosityA * uTangent * t_current * t_current' + viscosityB * uNormal * n_current * n_current' ) / rodParams.dt;
    
    Fd(index(1:2)) = Fd(index(1:2)) + fLocal/2;
    Jd(index(1:2),index(1:2)) = Jd(index(1:2),index(1:2)) + jac/2;
    
    
    velo = velo_2;
    
    uNormal = dot(n_current,velo);
    if (uNormal < 0)
        uNormal   = - uNormal;
        n_current = - n_current;
    end
    
    uTangent = dot(t_current,velo);
    if (uTangent < 0)
        uTangent  = - uTangent;
        t_current = - t_current;
    end
    
    fLocal = - 0.5 * 1000 * edge * (2 * rodParams.r0) * (viscosityA * uTangent * uTangent * t_current + viscosityB * uNormal * uNormal * n_current);
    
    jac = - 1000 * edge * (2 * rodParams.r0) * (viscosityA * uTangent * t_current * t_current' + viscosityB * uNormal * n_current * n_current' ) / rodParams.dt;
    
    Fd(index(3:4)) = Fd(index(3:4)) + fLocal/2;
    Jd(index(3:4),index(3:4)) = Jd(index(3:4),index(3:4)) + jac/2;
    
end

end
