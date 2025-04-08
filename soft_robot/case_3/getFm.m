function [Fm, Jm] = getFm(rodParams, sElement)
% This function computes the magnetic force and jacobian of the simulated system.
% Input:  rodParams - the defined beam struct contains the physical and
%                     numerical parameters of the simulated system
%
% Output: Fm - magnetic forces (ndof x 1)
%         Jm - magnetic jacobian (ndof x ndof)

Fm = zeros(rodParams.ndof, 1);
Jm = zeros(rodParams.ndof, rodParams.ndof);

for c=1:rodParams.ne
    
    brVec = sElement(c).Br;
    baVec = rodParams.Ba;
    
    node_1_start = sElement(c).start_node_1;
    node_2_start = sElement(c).start_node_2;
    
    t_start = (node_2_start - node_1_start) / norm(node_2_start - node_1_start);
    m_start = [-t_start(2); t_start(1)];
    
    node_1 = getVertex(rodParams.x0, sElement(c).nodeIndex(1));
    node_2 = getVertex(rodParams.x0, sElement(c).nodeIndex(2));
    
    t_current = (node_2 - node_1) / norm(node_2 - node_1);
    m_current = [-t_current(2); t_current(1)];
    
    edge = norm(node_2 - node_1);
    
    dtde = (eye(2) - t_current * t_current') / edge;
    dmde = -(t_current * m_current') / edge;
    
    dEde = sElement(c).refLen * rodParams.crossSection * ( dot(m_start, brVec) * dmde + dot(t_start, brVec) * dtde )' * baVec;
    
    dF = zeros(4,1);
    dF(1:2) = - dEde;
    dF(3:4) =   dEde;
    
    index = sElement(c).globalIndex;
    
    Fm(index) = Fm(index) - dF;
end

end
