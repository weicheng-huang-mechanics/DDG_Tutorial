function [dF, dJ] = stretchingForce(node0, node1, l_k, EA)

edge = (node1 - node0);
edgeLen = norm(edge);
tangent = edge / edgeLen;
epsX = edgeLen/l_k - 1;
dF_unit = EA * tangent * epsX;

dF = zeros(4,1);
dF(1:2) = - dF_unit;
dF(3:4) =   dF_unit;

Id3 = eye(2);
M = EA * ( (1/l_k - 1/edgeLen) * Id3 + 1/edgeLen * (edge*edge')/ edgeLen^2 );

dJ = zeros(4,4);
dJ(1:2, 1:2) = M;
dJ(3:4, 3:4) = M;
dJ(1:2, 3:4) = - M;
dJ(3:4, 1:2) = - M;

end
