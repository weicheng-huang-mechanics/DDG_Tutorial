% build stretching element
function sElement = InitialStretchingElement(rodParams, edge)

    EA = rodParams.EA;
    x = rodParams.x;
    ne = rodParams.ne;
    
    sElement = struct('nodeIndex', {}, 'globalIndex', {}, 'refLen', {}, 'EA_local', {}, ... 
        'start_node_1', {}, 'start_node_2', {}, 'Br', {});
    
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
        
        tangent = (node_2 - node_1) / norm(node_2 - node_1);
        
        sElement(i).Br = zeros(2,1);
        
        if ( abs(tangent(1)) < 1e-3 )
            
            xLocal = node_1(1);
            sElement(i).Br(1) = cos( 4 * pi * (xLocal-0.005)/0.07 );
            sElement(i).Br(2) = sin( 4 * pi * (xLocal-0.005)/0.07 );
        end    
    end
    
end