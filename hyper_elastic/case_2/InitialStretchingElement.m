function sElement = InitialStretchingElement(rodParams, edge)
% This function initializes stretching elements for a discrete rotational membrane.
%
%   Input:
%       rodParams - the defined rotational membrane struct contains the physical and
%                   numerical parameters of the simulated system
%       edge      - edge connectivity list (ne x 2)
%
%   Output:
%       sElement  - struct array containing stretching element information,
%                   including global indices, reference length, and stiffness

    E = rodParams.Y;
    nu = rodParams.nu;
    x = rodParams.x;
    ne = rodParams.ne;
    
    sElement = struct('nodeIndex', {}, 'globalIndex', {}, 'refLen', {}, 'refRadius', {}, 'refThickness', {}, 'E_local', {});
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
        
        sElement(i).refThickness = rodParams.r0;
        
        sElement(i).E_local = E / (2*(1+nu));
    end
    
end
