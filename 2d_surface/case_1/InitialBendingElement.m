function bElement = InitialBendingElement(rodParams, bend, sElement)

% This function initializes the bending elements for a rotational shell.
%
%   Input:
%       rodParams - the defined shell struct contains the physical and
%                   numerical parameters of the simulated system
%       bend      - bending element list (nb x 2)
%       sElement  - struct array representing all edges in the shell
%
%   Output:
%       bElement  - struct array of bending elements with 
%                   geometric and physical properties


EI = rodParams.EI;
x = rodParams.x;
nb = rodParams.nb;
nu = rodParams.nu;

bElement = struct('nodeIndex', {}, 'globalIndex', {}, 'nBar', {}, 'EI_local', {}, 'nu_local', {});

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
    
    xa = node_1(1);
    ya = node_1(2);
    xi = node_2(1);
    yi = node_2(2);
    xb = node_3(1);
    yb = node_3(2);
    l1 = sqrt((xi - xa)^2 + (yi - ya)^2);
    l2 = sqrt((xb - xi)^2 + (yb - yi)^2);

    t1x = (xi - xa)/l1;
    t1y = (yi - ya)/l1;
    t2x = (xb - xi)/l2;
    t2y = (yb - yi)/l2;
    tx = (t1x + t2x)/2;
    ty = (t1y + t2y)/2;

    bElement(i).kappaBar_1 = 2 *(t1x*t2y - t1y*t2x) / ((1 + t1x*t2x + t1y*t2y)*bElement(i).voroLen);
    bElement(i).kappaBar_2 = ty / (sqrt(tx*tx + ty*ty)*xi);
       
   
    % global index
    bElement(i).globalIndex = zeros(6,1);
    
    bElement(i).globalIndex(1) = 2 * (bElement(i).nodeIndex(1)-1) + 1;
    bElement(i).globalIndex(2) = 2 * (bElement(i).nodeIndex(1)-1) + 2;
    
    bElement(i).globalIndex(3) = 2 * (bElement(i).nodeIndex(2)-1) + 1;
    bElement(i).globalIndex(4) = 2 * (bElement(i).nodeIndex(2)-1) + 2;
    
    bElement(i).globalIndex(5) = 2 * (bElement(i).nodeIndex(3)-1) + 1;
    bElement(i).globalIndex(6) = 2 * (bElement(i).nodeIndex(3)-1) + 2;
    
    bElement(i).EI_local = EI;
    bElement(i).nu_local = nu;
    
end
