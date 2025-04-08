function [dF, dJ] = stretchingForce(node0, node1, l_k, EA)
% This function computes the stretching force and jacobian of a stretching element.
% Input:  node0 - position of the first node in the stretching element
%         node1 - position of the second node in the stretching element
%         l_k - reference length of the stretching element
%         EA - stretching stiffness
%
% Output: dF - stretching forces (6 x 1)
%         dJ - stretching jacobian (6 x 6)

edge = (node1 - node0);
edgeLen = norm(edge);
tangent = edge / edgeLen;
epsX = edgeLen/l_k - 1;
dF_unit = EA * tangent * epsX;

dF = zeros(6,1);
dF(1:3) = - dF_unit;
dF(4:6) =   dF_unit;

Id3 = eye(3);
M = EA * ( (1/l_k - 1/edgeLen) * Id3 + 1/edgeLen * (edge*edge')/ edgeLen^2 );

dJ = zeros(6,6);
dJ(1:3, 1:3) = M;
dJ(4:6, 4:6) = M;
dJ(1:3, 4:6) = - M;
dJ(4:6, 1:3) = - M;

end
