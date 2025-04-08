% build stretching element
function sElement = InitialStretchingElement(rodParams, edge)
% This function initializes stretching elements for a rotational shell.
%
%   Input:
%       rodParams - the defined shell struct contains the physical and
%                   numerical parameters of the simulated system
%       edge      - Edge connectivity list (ne x 2)
%
%   Output:
%       sElement  - Struct array containing stretching element information,
%                   including global indices, reference length, and stiffness

    EA = rodParams.EA;
    nu = rodParams.nu;
    x = rodParams.x;
    ne = rodParams.ne;
    
    sElement = struct('nodeIndex', {}, 'globalIndex', {}, 'refLen', {}, 'refRadius', {}, 'EA_local', {}, 'nu_local', {});
    for i = 1:ne
        sElement(i).nodeIndex = edge(i,:);
        
        % set global index
        sElement(i).globalIndex = zeros(4,1);
        
        sElement(i).globalIndex(1) = 2 * (sElement(i).nodeIndex(1) - 1) + 1;
        sElement(i).globalIndex(2) = 2 * (sElement(i).nodeIndex(1) - 1) + 2;
        
        sElement(i).globalIndex(3) = 2 * (sElement(i).nodeIndex(2) - 1) + 1;
        sElement(i).globalIndex(4) = 2 * (sElement(i).nodeIndex(2) - 1) + 2;
         
        node_1 = getVertex(x, sElement(i).nodeIndex(1));
        node_2 = getVertex(x, sElement(i).nodeIndex(2));
        
        sElement(i).refLen = norm(node_2 - node_1); 
        
        sElement(i).refRadius =  (node_2(1) + node_1(1)) / 2; 
        
        sElement(i).EA_local = EA;
        sElement(i).nu_local = nu;
    end
    
end
