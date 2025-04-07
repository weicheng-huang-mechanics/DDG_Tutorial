function sElement = InitialStretchingElement(plateParams, edge)
% This function initializes stretching elements for a discrete plate
%
%   Input:
%       plateParams - the defined plate struct contains the physical and
%                   numerical parameters of the simulated system
%       edge      - Edge connectivity list (ne x 2)
%
%   Output:
%       sElement  - Struct array containing stretching element information,
%                   including global indices, reference length, and stiffness

    EA = plateParams.EA;
    x  = plateParams.x;
    ne = plateParams.ne;
    
    sElement = struct('nodeIndex', {}, 'globalIndex', {}, 'refLen', {}, 'EA_local', {});

    for i = 1:ne
        sElement(i).nodeIndex = edge(i,:);
        
        % set global index
        sElement(i).globalIndex = zeros(6,1);
        
        sElement(i).globalIndex(1) = 3 * (sElement(i).nodeIndex(1) - 1) + 1;
        sElement(i).globalIndex(2) = 3 * (sElement(i).nodeIndex(1) - 1) + 2;
        sElement(i).globalIndex(3) = 3 * (sElement(i).nodeIndex(1) - 1) + 3;
        sElement(i).globalIndex(4) = 3 * (sElement(i).nodeIndex(2) - 1) + 1;
        sElement(i).globalIndex(5) = 3 * (sElement(i).nodeIndex(2) - 1) + 2;
        sElement(i).globalIndex(6) = 3 * (sElement(i).nodeIndex(2) - 1) + 3;

        node_1 = getVertex(x, sElement(i).nodeIndex(1));
        node_2 = getVertex(x, sElement(i).nodeIndex(2));
        
        sElement(i).refLen = norm(node_2 - node_1); 
        
        sElement(i).EA_local = EA;
    end
    
end