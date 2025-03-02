function [dF, dJ] = stretchingForce(node0, node1, l_k, r_k, EA)

edge = (node1 - node0);
edgeLen = norm(edge);
tangent = edge / edgeLen;
eps1 = edgeLen/l_k - 1;
dF_unit = EA * tangent * eps1;


radius = (node1(1) + node0(1)) / 2;

dF = zeros(4,1);

dF(1:2) = - dF_unit;
dF(3:4) =   dF_unit;

dF(1) = dF(1) + 0.5 * EA * (radius - r_k) / r_k;
dF(3) = dF(3) + 0.5 * EA * (radius - r_k) / r_k;

Id3 = eye(2);
M = EA * ( (1/l_k - 1/edgeLen) * Id3 + 1/edgeLen * (edge*edge')/ edgeLen^2 );

dJ = zeros(4,4);
dJ(1:2, 1:2) = M;
dJ(3:4, 3:4) = M;
dJ(1:2, 3:4) = - M;
dJ(3:4, 1:2) = - M;

dJ(1,1) = dJ(1,1) + 0.25 * EA / r_k;
dJ(1,3) = dJ(1,3) + 0.25 * EA / r_k;
dJ(3,1) = dJ(3,1) + 0.25 * EA / r_k;
dJ(3,3) = dJ(3,3) + 0.25 * EA / r_k;

dF = dF * 2 * pi * r_k;
dJ = dJ * 2 * pi * r_k;

end
