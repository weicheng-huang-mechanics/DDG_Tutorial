function [dF, dJ] = bendingForce(node0, node1, node2, kappaBar_1, kappaBar_2, refRadius, l_k, EI)

dF = zeros(6,1);
dJ = zeros(6,6);

m2e = [0;0;1];
m2f = [0;0;1];

gradKappa = zeros(6,1);

ee = (node1 - node0);
ef = (node2 - node1);

norm_e = norm(ee);
norm_f = norm(ef);

te = ee / norm_e;
tf = ef / norm_f;

kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;
tilde_d2 = (m2e + m2f) / chi;

kappa1 = 0.5 * dot( kb, m2e + m2f); 

Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + cross(tf,tilde_d2));
Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - cross(te,tilde_d2));

gradKappa(1:2) = -Dkappa1De(1:2);
gradKappa(3:4) = Dkappa1De(1:2) - Dkappa1Df(1:2);
gradKappa(5:6) = Dkappa1Df(1:2);


DDkappa1 = zeros(6, 6);

norm2_e = norm_e^2;
norm2_f = norm_f^2;

tt_o_tt = tilde_t' * tilde_t;
tmp = cross(tf, tilde_d2);
tf_c_d2t_o_tt = tmp' * tilde_t; 
tt_o_tf_c_d2t = tf_c_d2t_o_tt'; 
kb_o_d2e = kb' * m2e; 
d2e_o_kb = kb_o_d2e'; 

Id3 = eye(3);

D2kappa1De2 ...
    = 1.0 / norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t) ...
    - kappa1 / (chi * norm2_e) * (Id3 - te'*te) ...
    + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

tmp = cross(te, tilde_d2);
te_c_d2t_o_tt = tmp' * tilde_t;
tt_o_te_c_d2t = te_c_d2t_o_tt';
kb_o_d2f = kb' * m2f;
d2f_o_kb = kb_o_d2f';

D2kappa1Df2 ...
    = 1.0 / norm2_f * (2 * kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t) ...
    - kappa1 / (chi * norm2_f) * (Id3 - tf'*tf) ...
    + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

D2kappa1DeDf ...
    = -kappa1/(chi * norm_e * norm_f) * (Id3 + te'*tf) ...
    + 1.0 / (norm_e*norm_f) * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt + ...
    tt_o_te_c_d2t - crossMat(tilde_d2));
D2kappa1DfDe = D2kappa1DeDf';


DDkappa1(1:2, 1:2)  =   D2kappa1De2(1:2, 1:2);
DDkappa1(1:2, 3:4)  = - D2kappa1De2(1:2, 1:2) + D2kappa1DeDf(1:2, 1:2);
DDkappa1(1:2, 5:6)  =               - D2kappa1DeDf(1:2, 1:2);
DDkappa1(3:4, 1:2)  = - D2kappa1De2(1:2, 1:2)                + D2kappa1DfDe(1:2, 1:2);
DDkappa1(3:4, 3:4)  =   D2kappa1De2(1:2, 1:2) - D2kappa1DeDf(1:2, 1:2) - D2kappa1DfDe(1:2, 1:2) + D2kappa1Df2(1:2, 1:2);
DDkappa1(3:4, 5:6)  =                 D2kappa1DeDf(1:2, 1:2)                - D2kappa1Df2(1:2, 1:2);
DDkappa1(5:6, 1:2)  =                              - D2kappa1DfDe(1:2, 1:2);
DDkappa1(5:6, 3:4)  =                                D2kappa1DfDe(1:2, 1:2) - D2kappa1Df2(1:2, 1:2);
DDkappa1(5:6, 5:6)  =                                               D2kappa1Df2(1:2, 1:2);


dkappa = kappa1 - kappaBar_1;

dF = gradKappa * EI * dkappa / l_k;

dJ = 1.0 / l_k * gradKappa * EI * transpose(gradKappa);
temp = 1.0 / l_k * dkappa * EI;
dJ = dJ + temp * DDkappa1 ;


% grad and hession for kappa2

xa = node0(1);
ya = node0(2);
xi = node1(1);
yi = node1(2);
xb = node2(1);
yb = node2(2);

l = sqrt((xb - xa)^2 + (yb - ya)^2);
tx = (xb - xa) / l;
ty = (yb - ya) / l;
kappa2 = ty / (xi) - kappaBar_2;

dkappa2dx = zeros(6,1);

dkappa2dx(1) = ((2*xa - 2*xb)*(ya - yb))/(2*xi*((xa - xb)^2 + (ya - yb)^2)^(3/2));
dkappa2dx(2) = ((2*ya - 2*yb)*(ya - yb))/(2*xi*((xa - xb)^2 + (ya - yb)^2)^(3/2)) - 1/(xi*((xa - xb)^2 + (ya - yb)^2)^(1/2));
dkappa2dx(3) = (ya - yb)/(xi^2*((xa - xb)^2 + (ya - yb)^2)^(1/2));
dkappa2dx(4) = 0.0;
dkappa2dx(5) = -((2*xa - 2*xb)*(ya - yb))/(2*xi*((xa - xb)^2 + (ya - yb)^2)^(3/2));
dkappa2dx(6) = 1/(xi*((xa - xb)^2 + (ya - yb)^2)^(1/2)) - ((2*ya - 2*yb)*(ya - yb))/(2*xi*((xa - xb)^2 + (ya - yb)^2)^(3/2));

dkappa22dx2 = zeros(6,6);

dkappa22dx2(1,1) = ((ya - yb)*(- 2*xa^2 + 4*xa*xb - 2*xb^2 + ya^2 - 2*ya*yb + yb^2))/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));
dkappa22dx2(1,2) = -((xa - xb)*(- xa^2 + 2*xa*xb - xb^2 + 2*ya^2 - 4*ya*yb + 2*yb^2))/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));
dkappa22dx2(1,3) = -((2*xa - 2*xb)*(ya - yb))/(2*xi^2*((xa - xb)^2 + (ya - yb)^2)^(3/2));
dkappa22dx2(1,4) = 0.0;
dkappa22dx2(1,5) = -((ya - yb)*(- 2*xa^2 + 4*xa*xb - 2*xb^2 + ya^2 - 2*ya*yb + yb^2))/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));
dkappa22dx2(1,6) = ((xa - xb)*(- xa^2 + 2*xa*xb - xb^2 + 2*ya^2 - 4*ya*yb + 2*yb^2))/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));

dkappa22dx2(2,2) = ((xa - xb)*(- xa^2 + 2*xa*xb - xb^2 + 2*ya^2 - 4*ya*yb + 2*yb^2))/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));
dkappa22dx2(2,3) = ((2*xa - 2*xb)*(ya - yb))/(2*xi^2*((xa - xb)^2 + (ya - yb)^2)^(3/2));
dkappa22dx2(2,4) = 0.0;
dkappa22dx2(2,5) = ((ya - yb)*(- 2*xa^2 + 4*xa*xb - 2*xb^2 + ya^2 - 2*ya*yb + yb^2))/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));
dkappa22dx2(2,6) = -((xa - xb)*(- xa^2 + 2*xa*xb - xb^2 + 2*ya^2 - 4*ya*yb + 2*yb^2))/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));

dkappa22dx2(3,3) = -(2*(ya - yb))/(xi^3*((xa - xb)^2 + (ya - yb)^2)^(1/2));
dkappa22dx2(3,4) = 0.0;
dkappa22dx2(3,5) = ((2*xa - 2*xb)*(ya - yb))/(2*xi^2*((xa - xb)^2 + (ya - yb)^2)^(3/2));
dkappa22dx2(3,6) = -(xa - xb)^2/(xi^2*((xa - xb)^2 + (ya - yb)^2)^(3/2));

dkappa22dx2(4,4) = 0.0;
dkappa22dx2(4,5) = 0.0;
dkappa22dx2(4,6) = 0.0;

dkappa22dx2(5,5) = ((ya - yb)*(- 2*xa^2 + 4*xa*xb - 2*xb^2 + ya^2 - 2*ya*yb + yb^2))/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));
dkappa22dx2(5,6) = -((xa - xb)*(- xa^2 + 2*xa*xb - xb^2 + 2*ya^2 - 4*ya*yb + 2*yb^2))/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));

dkappa22dx2(6,6) = (3*ya*((xa - xb)^2 + (ya - yb)^2) - 3*yb*((xa - xb)^2 + (ya - yb)^2) - 9*ya*yb^2 + 9*ya^2*yb - 3*ya^3 + 3*yb^3)/(xi*((xa - xb)^2 + (ya - yb)^2)^(5/2));

dkappa22dx2(2,1) = dkappa22dx2(1,2);
dkappa22dx2(3,1) = dkappa22dx2(1,3);
dkappa22dx2(4,1) = dkappa22dx2(1,4);
dkappa22dx2(5,1) = dkappa22dx2(1,5);
dkappa22dx2(6,1) = dkappa22dx2(1,6);

dkappa22dx2(3,2) = dkappa22dx2(2,3);
dkappa22dx2(4,2) = dkappa22dx2(2,4);
dkappa22dx2(5,2) = dkappa22dx2(2,5);
dkappa22dx2(6,2) = dkappa22dx2(2,6);

dkappa22dx2(4,3) = dkappa22dx2(3,4);
dkappa22dx2(5,3) = dkappa22dx2(3,5);
dkappa22dx2(6,3) = dkappa22dx2(3,6);

dkappa22dx2(5,4) = dkappa22dx2(4,5);
dkappa22dx2(6,4) = dkappa22dx2(4,6);

dkappa22dx2(6,5) = dkappa22dx2(5,6);

dF = dF + kappa2 * dkappa2dx * EI * l_k;
dJ = dJ + (kappa2 * dkappa22dx2 + dkappa2dx * dkappa2dx') * EI * l_k ;

dF = dF * 2 * pi * refRadius;
dJ = dJ * 2 * pi * refRadius;
end
