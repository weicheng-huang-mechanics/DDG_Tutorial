function [Fp, Jp] = getFp(rodParams, sElement)
% This function computes the pressure force and jacobian of the simulated system.
% Input:  rodParams - the defined membrane struct contains the physical and
%                     numerical parameters of the simulated system
%         sElement - the stretching element list (ne x 2)
%
% Output: Fp - bending forces (2*nv x 1)
%         Jp - bending jacobian (2*nv x 2*nv)

Fp = zeros(rodParams.ndof, 1);
Jp = zeros(rodParams.ndof, rodParams.ndof);

totalF = zeros(4,1);

for c=1:rodParams.ne
    
    node_1 = getVertex(rodParams.x0, sElement(c).nodeIndex(1));
    node_2 = getVertex(rodParams.x0, sElement(c).nodeIndex(2));
    
    tangent = (node_2 - node_1) / norm(node_2 - node_1);
    normal = [-tangent(2);tangent(1)];
    
    edgeLen = norm(node_2 - node_1);
    
    pressureLoad = - normal * edgeLen * pi * (node_1(1) + node_2(1)) * rodParams.pressure / 2; 
    
    dF = zeros(4,1);
    
    dF(1:2) = pressureLoad;
    dF(3:4) = pressureLoad;
    
    totalF = totalF + dF;
    
    index = sElement(c).globalIndex;
    
    Fp(index) = Fp(index) - dF;
end

end
