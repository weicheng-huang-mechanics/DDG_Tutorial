function [dF, dJ] = stretchingForce(node0, node1, l0, r0, EA, nu)

xa = node0(1);
ya = node0(2);
xb = node1(1);
yb = node1(2);

ds = 2 * pi * l0 * r0;

l = sqrt( (xb-xa)*(xb-xa) + (yb-ya)*(yb-ya) );
r = (xa + xb) / 2;

eps1 = (l - l0) / l0;
eps2 = (r - r0) / r0;


deps1  = zeros(4,1);
deps2  = zeros(4,1);
ddeps1 = zeros(4,4);
ddeps2 = zeros(4,4);

deps1(1) = (2*xa - 2*xb)/(2*l0*((xa - xb)^2 + (ya - yb)^2)^(1/2));
deps1(2) = (2*ya - 2*yb)/(2*l0*((xa - xb)^2 + (ya - yb)^2)^(1/2));
deps1(3) = -(2*xa - 2*xb)/(2*l0*((xa - xb)^2 + (ya - yb)^2)^(1/2));
deps1(4) = -(2*ya - 2*yb)/(2*l0*((xa - xb)^2 + (ya - yb)^2)^(1/2));

deps2(1) = 1/(2*r0);
deps2(2) = 0;
deps2(3) = 1/(2*r0);
deps2(4) = 0;

ddeps1(1,1) = 1/(l0*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*xa - 2*xb)^2/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2));
ddeps1(1,2) = -((2*xa - 2*xb)*(2*ya - 2*yb))/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2));
ddeps1(1,3) = (2*xa - 2*xb)^2/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2)) - 1/(l0*((xa - xb)^2 + (ya - yb)^2)^(1/2));
ddeps1(1,4) = ((2*xa - 2*xb)*(2*ya - 2*yb))/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2));

ddeps1(2,2) = 1/(l0*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*ya - 2*yb)^2/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2));
ddeps1(2,3) = ((2*xa - 2*xb)*(2*ya - 2*yb))/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2));
ddeps1(2,4) = (2*ya - 2*yb)^2/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2)) - 1/(l0*((xa - xb)^2 + (ya - yb)^2)^(1/2));

ddeps1(3,3) = 1/(l0*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*xa - 2*xb)^2/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2));
ddeps1(3,4) = -((2*xa - 2*xb)*(2*ya - 2*yb))/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2));

ddeps1(4,4) = 1/(l0*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - (2*ya - 2*yb)^2/(4*l0*((xa - xb)^2 + (ya - yb)^2)^(3/2));

ddeps1(2,1) = ddeps1(1,2);
ddeps1(3,1) = ddeps1(1,3);
ddeps1(4,1) = ddeps1(1,4);
ddeps1(3,2) = ddeps1(2,3);
ddeps1(4,2) = ddeps1(2,4); 
ddeps1(4,3) = ddeps1(3,4);

dF = 2 * EA * ( eps1 * deps1 + eps2 * deps2 + nu * (eps1 * deps2 + eps2 * deps1) ) * ds;
dJ = 2 * EA * ( deps1 * deps1' + eps1 * ddeps1 + deps2 * deps2' + eps2 * ddeps2 + nu * (eps1 * ddeps2 + deps1 * deps2' + eps2 * ddeps1 + deps2 * deps1') ) * ds;

end
