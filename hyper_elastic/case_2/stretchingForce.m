function [dF, dJ] = stretchingForce(node0, node1, lBar, rBar, hBar, c1, c2)
% This function computes the stretching force and jacobian of a stretching element.
% Input:  node0 - position of the first node in the stretching element
%         node1 - position of the second node in the stretching element
%         l_k - reference length of the stretching element
%         EA - stretching stiffness
%
% Output: dF - stretching forces (4 x 1)
%         dJ - stretching jacobian (4 x 4)

dv = 2 * pi * lBar * rBar * hBar;

xa = node0(1);
ya = node0(2);
xb = node1(1);
yb = node1(2);

dF = zeros(4,1);
dJ = zeros(4,4);


dF(1) = c1*((2*xa - 2*xb)/lBar^2 + (xa/2 + xb/2)/rBar^2 - (lBar^2*rBar^2)/((xa/2 + xb/2)^3*((xa - xb)^2 + (ya - yb)^2)) - (lBar^2*rBar^2*(2*xa - 2*xb))/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2)) - c2*(rBar^2/(xa/2 + xb/2)^3 + (lBar^2*(2*xa - 2*xb))/((xa - xb)^2 + (ya - yb)^2)^2 - ((2*xa - 2*xb)*(xa/2 + xb/2)^2)/(lBar^2*rBar^2) - ((xa/2 + xb/2)*((xa - xb)^2 + (ya - yb)^2))/(lBar^2*rBar^2));
dF(2) = c1*((2*ya - 2*yb)/lBar^2 - (lBar^2*rBar^2*(2*ya - 2*yb))/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2)) - c2*((lBar^2*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^2 - ((xa/2 + xb/2)^2*(2*ya - 2*yb))/(lBar^2*rBar^2));
dF(3) = - c1*((2*xa - 2*xb)/lBar^2 - (xa/2 + xb/2)/rBar^2 + (lBar^2*rBar^2)/((xa/2 + xb/2)^3*((xa - xb)^2 + (ya - yb)^2)) - (lBar^2*rBar^2*(2*xa - 2*xb))/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2)) - c2*(rBar^2/(xa/2 + xb/2)^3 - (lBar^2*(2*xa - 2*xb))/((xa - xb)^2 + (ya - yb)^2)^2 + ((2*xa - 2*xb)*(xa/2 + xb/2)^2)/(lBar^2*rBar^2) - ((xa/2 + xb/2)*((xa - xb)^2 + (ya - yb)^2))/(lBar^2*rBar^2));
dF(4) = c2*((lBar^2*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^2 - ((xa/2 + xb/2)^2*(2*ya - 2*yb))/(lBar^2*rBar^2)) - c1*((2*ya - 2*yb)/lBar^2 - (lBar^2*rBar^2*(2*ya - 2*yb))/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2));

dJ(1,1) = c1*(2/lBar^2 + 1/(2*rBar^2) - (2*lBar^2*rBar^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2) + (3*lBar^2*rBar^2)/(2*(xa/2 + xb/2)^4*((xa - xb)^2 + (ya - yb)^2)) + (2*lBar^2*rBar^2*(2*xa - 2*xb))/((xa/2 + xb/2)^3*((xa - xb)^2 + (ya - yb)^2)^2) + (2*lBar^2*rBar^2*(2*xa - 2*xb)^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3)) + c2*((3*rBar^2)/(2*(xa/2 + xb/2)^4) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*(xa/2 + xb/2)^2)/(lBar^2*rBar^2) + (2*lBar^2*(2*xa - 2*xb)^2)/((xa - xb)^2 + (ya - yb)^2)^3 + ((xa - xb)^2 + (ya - yb)^2)/(2*lBar^2*rBar^2) + (2*(2*xa - 2*xb)*(xa/2 + xb/2))/(lBar^2*rBar^2));
dJ(1,2) = c2*(((xa/2 + xb/2)*(2*ya - 2*yb))/(lBar^2*rBar^2) + (2*lBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^3) + c1*((lBar^2*rBar^2*(2*ya - 2*yb))/((xa/2 + xb/2)^3*((xa - xb)^2 + (ya - yb)^2)^2) + (2*lBar^2*rBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3));
dJ(1,3) = c2*((2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (3*rBar^2)/(2*(xa/2 + xb/2)^4) - (2*(xa/2 + xb/2)^2)/(lBar^2*rBar^2) - (2*lBar^2*(2*xa - 2*xb)^2)/((xa - xb)^2 + (ya - yb)^2)^3 + ((xa - xb)^2 + (ya - yb)^2)/(2*lBar^2*rBar^2)) + c1*(1/(2*rBar^2) - 2/lBar^2 + (2*lBar^2*rBar^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2) + (3*lBar^2*rBar^2)/(2*(xa/2 + xb/2)^4*((xa - xb)^2 + (ya - yb)^2)) - (2*lBar^2*rBar^2*(2*xa - 2*xb)^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3));
dJ(1,4) = - c2*(((xa/2 + xb/2)*(2*ya - 2*yb))/(lBar^2*rBar^2) + (2*lBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^3) - c1*((lBar^2*rBar^2*(2*ya - 2*yb))/((xa/2 + xb/2)^3*((xa - xb)^2 + (ya - yb)^2)^2) + (2*lBar^2*rBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3));

dJ(2,2) = c2*((2*(xa/2 + xb/2)^2)/(lBar^2*rBar^2) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*lBar^2*(2*ya - 2*yb)^2)/((xa - xb)^2 + (ya - yb)^2)^3) + c1*(2/lBar^2 - (2*lBar^2*rBar^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2) + (2*lBar^2*rBar^2*(2*ya - 2*yb)^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3));
dJ(2,3) = c2*(((xa/2 + xb/2)*(2*ya - 2*yb))/(lBar^2*rBar^2) - (2*lBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^3) + c1*((lBar^2*rBar^2*(2*ya - 2*yb))/((xa/2 + xb/2)^3*((xa - xb)^2 + (ya - yb)^2)^2) - (2*lBar^2*rBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3));
dJ(2,4) = - c2*((2*(xa/2 + xb/2)^2)/(lBar^2*rBar^2) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*lBar^2*(2*ya - 2*yb)^2)/((xa - xb)^2 + (ya - yb)^2)^3) - c1*(2/lBar^2 - (2*lBar^2*rBar^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2) + (2*lBar^2*rBar^2*(2*ya - 2*yb)^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3));

dJ(3,3) = c1*(2/lBar^2 + 1/(2*rBar^2) - (2*lBar^2*rBar^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2) + (3*lBar^2*rBar^2)/(2*(xa/2 + xb/2)^4*((xa - xb)^2 + (ya - yb)^2)) - (2*lBar^2*rBar^2*(2*xa - 2*xb))/((xa/2 + xb/2)^3*((xa - xb)^2 + (ya - yb)^2)^2) + (2*lBar^2*rBar^2*(2*xa - 2*xb)^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3)) + c2*((3*rBar^2)/(2*(xa/2 + xb/2)^4) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*(xa/2 + xb/2)^2)/(lBar^2*rBar^2) + (2*lBar^2*(2*xa - 2*xb)^2)/((xa - xb)^2 + (ya - yb)^2)^3 + ((xa - xb)^2 + (ya - yb)^2)/(2*lBar^2*rBar^2) - (2*(2*xa - 2*xb)*(xa/2 + xb/2))/(lBar^2*rBar^2));
dJ(3,4) = - c2*(((xa/2 + xb/2)*(2*ya - 2*yb))/(lBar^2*rBar^2) - (2*lBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa - xb)^2 + (ya - yb)^2)^3) - c1*((lBar^2*rBar^2*(2*ya - 2*yb))/((xa/2 + xb/2)^3*((xa - xb)^2 + (ya - yb)^2)^2) - (2*lBar^2*rBar^2*(2*xa - 2*xb)*(2*ya - 2*yb))/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3));

dJ(4,4) = c2*((2*(xa/2 + xb/2)^2)/(lBar^2*rBar^2) - (2*lBar^2)/((xa - xb)^2 + (ya - yb)^2)^2 + (2*lBar^2*(2*ya - 2*yb)^2)/((xa - xb)^2 + (ya - yb)^2)^3) + c1*(2/lBar^2 - (2*lBar^2*rBar^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^2) + (2*lBar^2*rBar^2*(2*ya - 2*yb)^2)/((xa/2 + xb/2)^2*((xa - xb)^2 + (ya - yb)^2)^3));

dJ(2,1) = dJ(1,2);
dJ(3,1) = dJ(1,3);
dJ(4,1) = dJ(1,4);
dJ(3,2) = dJ(2,3);
dJ(4,2) = dJ(2,4);
dJ(4,3) = dJ(3,4);

dF = dF * dv;
dJ = dJ * dv;

end
