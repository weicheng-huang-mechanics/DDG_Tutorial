% build stretching element
function sElement = InitialStretchingElement(rodParams, edge)

EA = rodParams.EA;
x  = rodParams.x;
ne = rodParams.ne;

% build stretching element
sElement = struct('nodeIndex', {}, 'globalIndex', {}, 'refLen', {}, 't', {}, 'EA_local', {});

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
     
    nodePos_1 = getVertex(x, sElement(i).nodeIndex(1));
    nodePos_2 = getVertex(x, sElement(i).nodeIndex(2));
    
    % reference length
    sElement(i).refLen = norm(nodePos_2 - nodePos_1);
    
    % build tangent frame
    sElement(i).t = (nodePos_2 - nodePos_1) / norm(nodePos_2 - nodePos_1);
    
    sElement(i).EA_local = EA;
end
    
end