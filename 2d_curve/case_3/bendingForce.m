function [dF, dJ] = bendingForce(node0, node1, node2, m2e, m2f, ...
    kappaBar, l_k, EI)
% This function computes the bending force and jacobian of a bending
% element.
% INPUTS: node0 - position of the first node in the bending element
%         node1 - position of the second node in the bending element
%         node2 - poisiton of the third node in the bending element
%         m2e - material frame of the first edge 
%         m2f - material frame of the second edge
%         kappaBar - the natural curvature of the bending element
%         l_k - voronoi length of the bending element
%         EI - bending stiffness
%
% OUTPUTS: dF - bending forces (6 x 1)
%          dJ - bending jacobian (6 x 6)

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


dkappa = kappa1 - kappaBar;

dF = gradKappa * EI * dkappa / l_k;

dJ = 1.0 / l_k * gradKappa * EI * transpose(gradKappa);
temp = 1.0 / l_k * dkappa * EI;
dJ = dJ + temp * DDkappa1 ;

end
