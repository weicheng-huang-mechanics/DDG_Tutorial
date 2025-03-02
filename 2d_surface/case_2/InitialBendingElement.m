function bElement = InitialBendingElement(rodParams, bend, sElement)

EI = rodParams.EI;
x = rodParams.x;
nb = rodParams.nb;

bElement = struct('nodeIndex', {}, 'globalIndex', {}, 'nBar', {}, 'EI_local', {});

for i = 1:nb
    bElement(i).edgeIndex = bend(i,:);  
    
    % edge used in bending
    localEdge_1 = sElement(bElement(i).edgeIndex(1));
    localEdge_2 = sElement(bElement(i).edgeIndex(2));   
    
    % local node index
    index1 = localEdge_1.nodeIndex(1);
    index2 = localEdge_1.nodeIndex(2);
    index3 = localEdge_2.nodeIndex(1);
    index4 = localEdge_2.nodeIndex(2);
    
    if (index1 == index3)
        bElement(i).nodeIndex(1) = index2;
        bElement(i).nodeIndex(2) = index1;
        bElement(i).nodeIndex(3) = index4;
    end
    
    if (index1 == index4)
        bElement(i).nodeIndex(1) = index2;
        bElement(i).nodeIndex(2) = index1;
        bElement(i).nodeIndex(3) = index3;
    end
    
    if (index2 == index3)
        bElement(i).nodeIndex(1) = index1;
        bElement(i).nodeIndex(2) = index2;
        bElement(i).nodeIndex(3) = index4;
    end
    
    if (index2 == index4)
        bElement(i).nodeIndex(1) = index1;
        bElement(i).nodeIndex(2) = index2;
        bElement(i).nodeIndex(3) = index3;
    end
    
    node_1 = getVertex(x, bElement(i).nodeIndex(1));
    node_2 = getVertex(x, bElement(i).nodeIndex(2));
    node_3 = getVertex(x, bElement(i).nodeIndex(3));
    
    % 
    bElement(i).voroLen = ( norm(node_3 - node_2) + norm(node_2 - node_1) ) / 2;
    bElement(i).refRadius = node_2(1);
    
    te(1:2,1) = (node_2 - node_1) / norm(node_2 - node_1);
    tf(1:2,1) = (node_3 - node_2) / norm(node_3 - node_2);
    te(3,1)   = 0.0;
    tf(3,1)   = 0.0;
    
    kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));
    
    % kappaBar
    bElement(i).kappaBar_1 = dot(kb, [0;0;1] ); 
    
    l = norm(node_3 - node_1);
    ty = (node_3(2) - node_1(2)) / l;

    bElement(i).kappaBar_2 = ty / node_2(1);
   
   
    % global index
    bElement(i).globalIndex = zeros(6,1);
    
    bElement(i).globalIndex(1) = 2 * (bElement(i).nodeIndex(1)-1) + 1;
    bElement(i).globalIndex(2) = 2 * (bElement(i).nodeIndex(1)-1) + 2;
    
    bElement(i).globalIndex(3) = 2 * (bElement(i).nodeIndex(2)-1) + 1;
    bElement(i).globalIndex(4) = 2 * (bElement(i).nodeIndex(2)-1) + 2;
    
    bElement(i).globalIndex(5) = 2 * (bElement(i).nodeIndex(3)-1) + 1;
    bElement(i).globalIndex(6) = 2 * (bElement(i).nodeIndex(3)-1) + 2;
    
    bElement(i).EI_local = EI;
    
end