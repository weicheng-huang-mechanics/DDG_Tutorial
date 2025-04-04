% build stretching element
function sElement = InitialStretchingElement(rodParams, edge)

    EA = rodParams.EA;
    x = rodParams.x;
    ne = rodParams.ne;
    
    sElement = struct('nodeIndex', {}, 'globalIndex', {}, 'refLen', {}, 'EA_local', {}, 'c1', {}, 'c2', {});
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
        
        sElement(i).EA_local = EA;
        
        sElement(i).c1 = 0.8 * EA / 6;
        sElement(i).c2 = 0.2 * EA / 6;
    end
    
end