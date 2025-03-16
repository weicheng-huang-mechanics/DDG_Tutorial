% build stretching element
function sElement = InitialStretchingElement(rodParams, edge)

    EA = rodParams.EA;
    x = rodParams.x;
    ne = rodParams.ne;
    
    sElement = struct('nodeIndex', {}, 'globalIndex', {}, 'refLen', {}, 'EA_local', {}, ... 
        'start_node_1', {}, 'start_node_2', {}, 'Br', {});
    
    figure(3);
    hold on
    
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
        
        sElement(i).start_node_1 = node_1;
        sElement(i).start_node_2 = node_2;
        
        sElement(i).Br = zeros(2,1);
        
        sElement(i).Br = [1e1; 0.0];    
        
        sElement(i).Br(1) = 1e2 * cos(2 * pi * (i-1) / (ne-1) );
        sElement(i).Br(2) = 1e2 * sin(2 * pi * (i-1) / (ne-1) );
    end
    
    
end
