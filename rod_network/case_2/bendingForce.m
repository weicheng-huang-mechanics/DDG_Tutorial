function [dF, dJ] = bendingForce(node0, node1, node2, m1e, m2e, m1f, m2f, ...
    kappaBar, l_k, EI, sign_1, sign_2)
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
% OUTPUTS: dF - bending forces (11 x 1)
%          dJ - bending jacobian (11 x 11)

gradKappa = zeros(11,2);

ee = (node1 - node0);
ef = (node2 - node1);

norm_e = norm(ee);
norm_f = norm(ef);

te = ee / norm_e;
tf = ef / norm_f;

kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;
tilde_d1 = (m1e + m1f) / chi;
tilde_d2 = (m2e + m2f) / chi;


kappa1 = 0.5 * dot( kb, m2e + m2f); 
kappa2 = -0.5 * dot( kb, m1e + m1f); 

Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + cross(tf,tilde_d2));
Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - cross(te,tilde_d2));

Dkappa2De = 1.0 / norm_e * (-kappa2 * tilde_t - cross(tf,tilde_d1));
Dkappa2Df = 1.0 / norm_f * (-kappa2 * tilde_t + cross(te,tilde_d1));

gradKappa(1:3, 1) = -Dkappa1De;
gradKappa(5:7, 1) = Dkappa1De - Dkappa1Df;
gradKappa(9:11, 1) = Dkappa1Df;

gradKappa(1:3, 2) = -Dkappa2De;
gradKappa(5:7, 2) = Dkappa2De - Dkappa2Df;
gradKappa(9:11, 2) = Dkappa2Df;

gradKappa(4, 1) = -0.5 * dot(kb, m1e);
gradKappa(8, 1) = -0.5 * dot(kb, m1f);
gradKappa(4, 2) = -0.5 * dot(kb, m2e);
gradKappa(8, 2) = -0.5 * dot(kb, m2f);

if (sign_1 < 0)
    gradKappa(4, 1) = - gradKappa(4, 1);
    gradKappa(4, 2) = - gradKappa(4, 2);
end

if (sign_2 < 0)
    gradKappa(8, 1) = - gradKappa(8, 1);
    gradKappa(8, 2) = - gradKappa(8, 2);
end

DDkappa1 = zeros(11, 11);
DDkappa2 = zeros(11, 11);

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

tmp = cross(tf, tilde_d1);
tf_c_d1t_o_tt = tmp'*tilde_t; 
tt_o_tf_c_d1t = tf_c_d1t_o_tt'; 
kb_o_d1e = kb'*m1e; 
d1e_o_kb = kb_o_d1e'; 

D2kappa2De2 ...
    = 1.0 / norm2_e * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt + tt_o_tf_c_d1t) ...
    - kappa2 / (chi * norm2_e) * (Id3 - te'*te) ...
    - 1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);

tmp = cross(te, tilde_d1);
te_c_d1t_o_tt = tmp'*tilde_t; 
tt_o_te_c_d1t = te_c_d1t_o_tt'; 
kb_o_d1f = kb'*m1f; 
d1f_o_kb =  kb_o_d1f'; 

D2kappa2Df2 ...
    = 1.0 / norm2_f * (2 * kappa2 * tt_o_tt - te_c_d1t_o_tt - tt_o_te_c_d1t) ...
    - kappa2 / (chi * norm2_f) * (Id3 - tf'*tf) ...
    - 1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb);

D2kappa2DeDf ...
    = -kappa2/(chi * norm_e * norm_f) * (Id3 + te'*tf) ...
    + 1.0 / (norm_e*norm_f) * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t + crossMat(tilde_d1));

D2kappa2DfDe = D2kappa2DeDf';

D2kappa1Dthetae2 = -0.5 * dot(kb, m2e);
D2kappa1Dthetaf2 = -0.5 * dot(kb, m2f);
D2kappa2Dthetae2 =  0.5 * dot(kb, m1e);
D2kappa2Dthetaf2 =  0.5 * dot(kb, m1f);

D2kappa1DeDthetae ...
    = 1.0 / norm_e * (0.5 * dot(kb, m1e) * tilde_t - 1.0 / chi * cross(tf, m1e));
D2kappa1DeDthetaf ...
    = 1.0 / norm_e * (0.5 * dot(kb, m1f) * tilde_t - 1.0 / chi * cross(tf, m1f));
D2kappa1DfDthetae ...
    = 1.0 / norm_f * (0.5 * dot(kb, m1e) * tilde_t + 1.0 / chi * cross(te, m1e));
D2kappa1DfDthetaf ...
    = 1.0 / norm_f * (0.5 * dot(kb, m1f) * tilde_t + 1.0 / chi * cross(te, m1f));
D2kappa2DeDthetae ...
    = 1.0 / norm_e * (0.5 * dot(kb, m2e) * tilde_t - 1.0 / chi * cross(tf, m2e));
D2kappa2DeDthetaf ...
    = 1.0 / norm_e * (0.5 * dot(kb, m2f) * tilde_t - 1.0 / chi * cross(tf, m2f));
D2kappa2DfDthetae ...
    = 1.0 / norm_f * (0.5 * dot(kb, m2e) * tilde_t + 1.0 / chi * cross(te, m2e));
D2kappa2DfDthetaf ...
    = 1.0 / norm_f * (0.5 * dot(kb, m2f) * tilde_t + 1.0 / chi * cross(te, m2f));

DDkappa1(1:3, 1:3)  =   D2kappa1De2;
DDkappa1(1:3, 5:7)  = - D2kappa1De2 + D2kappa1DeDf;
DDkappa1(1:3, 9:11) =               - D2kappa1DeDf;
DDkappa1(5:7, 1:3)  = - D2kappa1De2                + D2kappa1DfDe;
DDkappa1(5:7, 5:7)  =   D2kappa1De2 - D2kappa1DeDf - D2kappa1DfDe + D2kappa1Df2;
DDkappa1(5:7, 9:11) =                 D2kappa1DeDf                - D2kappa1Df2;
DDkappa1(9:11, 1:3)  =                              - D2kappa1DfDe;
DDkappa1(9:11, 5:7)  =                                D2kappa1DfDe - D2kappa1Df2;
DDkappa1(9:11, 9:11) =                                               D2kappa1Df2;

DDkappa1(4, 4)     =   D2kappa1Dthetae2;
DDkappa1(8, 8)     =   D2kappa1Dthetaf2;

DDkappa1(1:3, 4)   = - D2kappa1DeDthetae;
DDkappa1(5:7, 4)   =   D2kappa1DeDthetae - D2kappa1DfDthetae;
DDkappa1(9:11,4)   =                       D2kappa1DfDthetae;
DDkappa1(4, 1:3)   =   transpose(DDkappa1(1:3, 4));
DDkappa1(4, 5:7)   =   transpose(DDkappa1(5:7, 4));
DDkappa1(4, 9:11)  =   transpose(DDkappa1(9:11,4));

DDkappa1(1:3, 8)   = - D2kappa1DeDthetaf;
DDkappa1(5:7, 8)   =   D2kappa1DeDthetaf - D2kappa1DfDthetaf;
DDkappa1(9:11, 8)  =                       D2kappa1DfDthetaf;
DDkappa1(8, 1:3)   =   transpose(DDkappa1(1:3, 8));
DDkappa1(8, 5:7)   =   transpose(DDkappa1(5:7, 8));
DDkappa1(8, 9:11)  =   transpose(DDkappa1(9:11,8));

DDkappa2(1:3, 1:3) =   D2kappa2De2;
DDkappa2(1:3, 5:7) = - D2kappa2De2 + D2kappa2DeDf;
DDkappa2(1:3, 9:11) =               - D2kappa2DeDf;
DDkappa2(5:7, 1:3) = - D2kappa2De2                + D2kappa2DfDe;
DDkappa2(5:7, 5:7) =   D2kappa2De2 - D2kappa2DeDf - D2kappa2DfDe + D2kappa2Df2;
DDkappa2(5:7, 9:11)=                 D2kappa2DeDf                - D2kappa2Df2;
DDkappa2(9:11, 1:3)=                              - D2kappa2DfDe;
DDkappa2(9:11, 5:7)=                                D2kappa2DfDe - D2kappa2Df2;
DDkappa2(9:11, 9:11)=                                               D2kappa2Df2;

DDkappa2(4, 4)     = D2kappa2Dthetae2;
DDkappa2(8, 8)     = D2kappa2Dthetaf2;

DDkappa2(1:3, 4)   = - D2kappa2DeDthetae;
DDkappa2(5:7, 4)   =   D2kappa2DeDthetae - D2kappa2DfDthetae;
DDkappa2(9:11,4)   =                       D2kappa2DfDthetae;
DDkappa2(4, 1:3)   =   transpose(DDkappa2(1:3, 4));
DDkappa2(4, 5:7)   =   transpose(DDkappa2(5:7, 4));
DDkappa2(4, 9:11)  =   transpose(DDkappa2(9:11,4));

DDkappa2(1:3, 8)   = - D2kappa2DeDthetaf;
DDkappa2(5:7, 8)   =   D2kappa2DeDthetaf - D2kappa2DfDthetaf;
DDkappa2(9:11,8)   =                       D2kappa2DfDthetaf;
DDkappa2(8, 1:3)   =   transpose(DDkappa2(1:3, 8));
DDkappa2(8, 5:7)   =   transpose(DDkappa2(5:7, 8));
DDkappa2(8,9:11)   =   transpose(DDkappa2(9:11,8));

if (sign_1 < 0)
    DDkappa1(1:3,  4)  = - DDkappa1(1:3,  4);
    DDkappa1(5:7,  4)  = - DDkappa1(5:7,  4);
    DDkappa1(9:11, 4)  = - DDkappa1(9:11, 4);
    DDkappa1(1:3,  8)  = - DDkappa1(1:3,  8);
    DDkappa1(5:7,  8)  = - DDkappa1(5:7,  8);
    DDkappa1(9:11, 8)  = - DDkappa1(9:11, 8);
    
    DDkappa1(4, 1:3)   = - DDkappa1(4, 1:3);
    DDkappa1(4, 5:7)   = - DDkappa1(4, 5:7);
    DDkappa1(4, 9:11)  = - DDkappa1(4, 9:11);
    DDkappa1(8, 1:3)   = - DDkappa1(8, 1:3);
    DDkappa1(8, 5:7)   = - DDkappa1(8, 5:7);
    DDkappa1(8, 9:11)  = - DDkappa1(8, 9:11);
end

if (sign_2 < 0)
    DDkappa2(1:3,  4)  = - DDkappa2(1:3,  4);
    DDkappa2(5:7,  4)  = - DDkappa2(5:7,  4);
    DDkappa2(9:11, 4)  = - DDkappa2(9:11, 4);
    DDkappa2(1:3,  8)  = - DDkappa2(1:3,  8);
    DDkappa2(5:7,  8)  = - DDkappa2(5:7,  8);
    DDkappa2(9:11, 8)  = - DDkappa2(9:11, 8);
    
    DDkappa2(4, 1:3)   = - DDkappa2(4, 1:3);
    DDkappa2(4, 5:7)   = - DDkappa2(4, 5:7);
    DDkappa2(4, 9:11)  = - DDkappa2(4, 9:11);
    DDkappa2(8, 1:3)   = - DDkappa2(8, 1:3);
    DDkappa2(8, 5:7)   = - DDkappa2(8, 5:7);
    DDkappa2(8, 9:11)  = - DDkappa2(8, 9:11);
end

EIMat = [EI 0; ...
    0 EI];
kappaVector = [kappa1 kappa2];

dkappaVector = kappaVector - kappaBar;
dF = gradKappa * EIMat * dkappaVector' / l_k;

dJ = 1.0 / l_k * gradKappa * EIMat * transpose(gradKappa);
temp = 1.0 / l_k * dkappaVector * EIMat;
dJ = dJ + temp(1) * DDkappa1 + temp(2) * DDkappa2;

end
