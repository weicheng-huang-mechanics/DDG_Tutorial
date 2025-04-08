function [dF, dJ] = stretchingForce(node0, node1, lBar, c1, c2)
% This function computes the stretching force and jacobian of a stretching element.
% Input:  node0 - position of the first node in the stretching element
%         node1 - position of the second node in the stretching element
%         l_k - voronoi length of the bending element
%         EA - stretching stiffness
%
% Output: dF - stretching forces (4 x 1)
%         dJ - stretching jacobian (4 x 4)

xa = node0(1);
ya = node0(2);
xb = node1(1);
yb = node1(2);

dF = zeros(4,1);
dJ = zeros(4,4);

dF(1) = c1*((2*xa - 2*xb)/lBar^2 - (lBar*(2*xa - 2*xb))/((xa - xb)^2 + (ya - yb)^2)^(3/2)) + c2*((2*xa - 2*xb)/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (lBar^2*(2*xa - 2*xb))/((xa - xb)^2 + (ya - yb)^2)^2);
dF(2) = c1*((2*ya - 2*yb)/lBar^2 - (lBar*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^(3/2)) + c2*((2*ya - 2*yb)/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (lBar^2*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^2);
dF(3) = - c1*((2*xa - 2*xb)/lBar^2 - (lBar*(2*xa - 2*xb))/((xa - xb)^2 + (ya - yb)^2)^(3/2)) - c2*((2*xa - 2*xb)/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (lBar^2*(2*xa - 2*xb))/((xa - xb)^2 + (ya - yb)^2)^2);
dF(4) = - c1*((2*ya - 2*yb)/lBar^2 - (lBar*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^(3/2)) - c2*((2*ya - 2*yb)/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (lBar^2*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^2);

dJ(1,1) = c1*(2/lBar^2 - (2*lBar)/((xa - xb)^2 + (ya - yb)^2)^(3/2) + (3*lBar*(2*xa - 2*xb)^2)/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2))) + c2*(2/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*lBar^2*(2*xa - 2*xb)^2)/((xa - xb)^2 + (ya - yb)^2)^3 - (2*xa - 2*xb)^2/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2)));
dJ(1,2) = c2*((2*lBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^3 - ((2*xa - 2*xb)*(2*ya - 2*yb))/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2))) + (3*c1*lBar*(2*xa - 2*xb)*(2*ya - 2*yb))/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2));
dJ(1,3) = - c1*(2/lBar^2 - (2*lBar)/((xa - xb)^2 + (ya - yb)^2)^(3/2) + (3*lBar*(2*xa - 2*xb)^2)/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2))) - c2*(2/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*lBar^2*(2*xa - 2*xb)^2)/((xa - xb)^2 + (ya - yb)^2)^3 - (2*xa - 2*xb)^2/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2)));
dJ(1,4) = - c2*((2*lBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^3 - ((2*xa - 2*xb)*(2*ya - 2*yb))/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2))) - (3*c1*lBar*(2*xa - 2*xb)*(2*ya - 2*yb))/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2));

dJ(2,2) = c1*(2/lBar^2 - (2*lBar)/((xa - xb)^2 + (ya - yb)^2)^(3/2) + (3*lBar*(2*ya - 2*yb)^2)/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2))) + c2*(2/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*lBar^2*(2*ya - 2*yb)^2)/((xa - xb)^2 + (ya - yb)^2)^3 - (2*ya - 2*yb)^2/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2)));
dJ(2,3) = - c2*((2*lBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^3 - ((2*xa - 2*xb)*(2*ya - 2*yb))/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2))) - (3*c1*lBar*(2*xa - 2*xb)*(2*ya - 2*yb))/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2));
dJ(2,4) = - c1*(2/lBar^2 - (2*lBar)/((xa - xb)^2 + (ya - yb)^2)^(3/2) + (3*lBar*(2*ya - 2*yb)^2)/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2))) - c2*(2/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*lBar^2*(2*ya - 2*yb)^2)/((xa - xb)^2 + (ya - yb)^2)^3 - (2*ya - 2*yb)^2/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2)));
 
dJ(3,3) = c1*(2/lBar^2 - (2*lBar)/((xa - xb)^2 + (ya - yb)^2)^(3/2) + (3*lBar*(2*xa - 2*xb)^2)/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2))) + c2*(2/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*lBar^2*(2*xa - 2*xb)^2)/((xa - xb)^2 + (ya - yb)^2)^3 - (2*xa - 2*xb)^2/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2)));
dJ(3,4) = c2*((2*lBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^3 - ((2*xa - 2*xb)*(2*ya - 2*yb))/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2))) + (3*c1*lBar*(2*xa - 2*xb)*(2*ya - 2*yb))/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2));

dJ(4,4) = c1*(2/lBar^2 - (2*lBar)/((xa - xb)^2 + (ya - yb)^2)^(3/2) + (3*lBar*(2*ya - 2*yb)^2)/(2*((xa - xb)^2 + (ya - yb)^2)^(5/2))) + c2*(2/(lBar*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*lBar^2*(2*ya - 2*yb)^2)/((xa - xb)^2 + (ya - yb)^2)^3 - (2*ya - 2*yb)^2/(2*lBar*((xa - xb)^2 + (ya - yb)^2)^(3/2)));

dJ(2,1) = dJ(1,2);
dJ(3,1) = dJ(1,3);
dJ(4,1) = dJ(1,4);
dJ(3,2) = dJ(2,3);
dJ(4,2) = dJ(2,4);
dJ(4,3) = dJ(3,4);

dF = dF * lBar;
dJ = dJ * lBar;

end
