function [Fp, Jp] = getFp(rodParams, sElement)

% This function computes the pressure force of the simulated system.
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%
% OUTPUTS: Fp - gravitational forces (2*nv x 1)
%          Jp - gravitational jacobian (2*nv x 2*nv)

Fp = zeros(rodParams.ndof, 1);
Jp = zeros(rodParams.ndof, rodParams.ndof);

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
    
    index = sElement(c).globalIndex;
    
    Fp(index) = Fp(index) - dF;
end


end
