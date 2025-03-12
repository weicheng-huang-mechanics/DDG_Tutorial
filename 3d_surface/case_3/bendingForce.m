function [dF, dJ] = bendingForce(x_1, x_2, x_3, x_4, nBar_Scale, EI)

% Compute Edge
e_1 = x_3 - x_1;
e_2 = x_2 - x_1;
e_3 = x_2 - x_4;
e_4 = x_3 - x_4;

% Compute normal vectors
n_1 = cross(e_1, e_2);
n_2 = cross(e_3, e_4);

% Normalize the normal vectors
norm_1 = n_1 / norm(n_1);
norm_2 = n_2 / norm(n_2);

if ( norm((norm_1 - norm_2) ) > 1e-9)
    norm_12 = (norm_1 - norm_2) / norm((norm_1 - norm_2) );
else
    norm_12 = (norm_1 - norm_2) / (norm((norm_1 - norm_2) )+1e-9);
end

nCurrent = norm(norm_1 - norm_2) - nBar_Scale;

% Identity matrix
Id3 = eye(3);

% Compute gradient matrices
gradN1 = (Id3 - (norm_1 * norm_1')) / norm(n_1);
gradN2 = (Id3 - (norm_2 * norm_2')) / norm(n_2);

% Compute force contributions
% dEde1 =   (norm_1 - norm_2)' * gradN1 * crossMat(e_2);
% dEde2 = - (norm_1 - norm_2)' * gradN1 * crossMat(e_1);
% dEde3 = - (norm_1 - norm_2)' * gradN2 * crossMat(e_4);
% dEde4 =   (norm_1 - norm_2)' * gradN2 * crossMat(e_3);

dEde1 =   nCurrent * norm_12' * gradN1 * crossMat(e_2);
dEde2 = - nCurrent * norm_12' * gradN1 * crossMat(e_1);
dEde3 = - nCurrent * norm_12' * gradN2 * crossMat(e_4);
dEde4 =   nCurrent * norm_12' * gradN2 * crossMat(e_3);

% Initialize force vector
dF = zeros(12,1);

% Assign force components
dF(1:3)   = - dEde1 - dEde2;
dF(4:6)   =   dEde3 + dEde2;
dF(7:9)   =   dEde1 + dEde4;
dF(10:12) = - dEde4 - dEde3;

dF = - dF * EI;

% Compute second derivatives (Jacobian components)
d2Ede12 = -crossMat(e_2) * gradN1 * gradN1 * crossMat(e_2) ...
    - crossMat(e_2) * ( ( nCurrent * norm_12(1) ) * hessionMatrix_1(n_1) ...
    + ( nCurrent * norm_12(2) ) * hessionMatrix_2(n_1) ...
    + ( nCurrent * norm_12(3) ) * hessionMatrix_3(n_1) ) * crossMat(e_2);

d2Ede22 = -crossMat(e_1) * gradN1 * gradN1 * crossMat(e_1) ...
    - crossMat(e_1) * ( ( nCurrent * norm_12(1) ) * hessionMatrix_1(n_1) ...
    + ( nCurrent * norm_12(2) ) * hessionMatrix_2(n_1) ...
    + ( nCurrent * norm_12(3) ) * hessionMatrix_3(n_1) ) * crossMat(e_1);

d2Ede32 = -crossMat(e_4) * gradN2 * gradN2 * crossMat(e_4) ...
    + crossMat(e_4) * ( ( nCurrent * norm_12(1) ) * hessionMatrix_1(n_2) ...
    + ( nCurrent * norm_12(2) ) * hessionMatrix_2(n_2) ...
    + ( nCurrent * norm_12(3) ) * hessionMatrix_3(n_2) ) * crossMat(e_4);

d2Ede42 = -crossMat(e_3) * gradN2 * gradN2 * crossMat(e_3) ...
    + crossMat(e_3) * ( ( nCurrent * norm_12(1) ) * hessionMatrix_1(n_2) ...
    + ( nCurrent * norm_12(2) ) * hessionMatrix_2(n_2) ...
    + ( nCurrent * norm_12(3) ) * hessionMatrix_3(n_2) ) * crossMat(e_3);

% Compute mixed derivatives
d2Ede1de2 = crossMat(e_1) * gradN1 * gradN1 * crossMat(e_2) ...
    + crossMat(e_1) * ( ( nCurrent * norm_12(1) ) * hessionMatrix_1(n_1) ...
    + ( nCurrent * norm_12(2) ) * hessionMatrix_2(n_1) ...
    + ( nCurrent * norm_12(3) ) * hessionMatrix_3(n_1) ) * crossMat(e_2);
d2Ede2de1 = d2Ede1de2';

d2Ede1de3 = crossMat(e_2) * gradN1 * gradN2 * crossMat(e_4);
d2Ede3de1 = d2Ede1de3';

d2Ede1de4 = -crossMat(e_2) * gradN1 * gradN2 * crossMat(e_3);
d2Ede4de1 = d2Ede1de4';

d2Ede2de3 = -crossMat(e_1) * gradN1 * gradN2 * crossMat(e_4);
d2Ede3de2 = d2Ede2de3';

d2Ede2de4 = crossMat(e_1) * gradN1 * gradN2 * crossMat(e_3);
d2Ede4de2 = d2Ede2de4';

d2Ede3de4 = crossMat(e_3) * gradN2 * gradN2 * crossMat(e_4) ...
    - crossMat(e_3) * ( ( nCurrent * norm_12(1) ) * hessionMatrix_1(n_2) ...
    + ( nCurrent * norm_12(2) ) * hessionMatrix_2(n_2) ...
    + ( nCurrent * norm_12(3) ) * hessionMatrix_3(n_2) ) * crossMat(e_4);
d2Ede4de3 = d2Ede3de4';

% Initialize Jacobian matrix
dJ = zeros(12,12);

% Assign Jacobian blocks
dJ(1:3,1:3)     = d2Ede12 + d2Ede22 + d2Ede1de2 + d2Ede2de1;
dJ(4:6,4:6)     = d2Ede22 + d2Ede32 + d2Ede2de3 + d2Ede3de2;
dJ(7:9,7:9)     = d2Ede12 + d2Ede42 + d2Ede1de4 + d2Ede4de1;
dJ(10:12,10:12) = d2Ede32 + d2Ede42 + d2Ede3de4 + d2Ede4de3;

dJ(1:3,4:6) = - d2Ede22 - d2Ede1de2 - d2Ede2de3 - d2Ede1de3;
dJ(4:6,1:3) = dJ(1:3,4:6)';

dJ(1:3,7:9) = - d2Ede12 - d2Ede1de4 - d2Ede2de1 - d2Ede2de4;
dJ(7:9,1:3) = dJ(1:3,7:9)';

dJ(1:3,10:12) = d2Ede1de3 + d2Ede1de4 + d2Ede2de3 + d2Ede2de4;
dJ(10:12,1:3) = dJ(1:3,10:12)';

dJ(4:6,7:9) = d2Ede2de1 + d2Ede2de4 + d2Ede3de1 + d2Ede3de4;
dJ(7:9,4:6) = dJ(4:6,7:9)';

dJ(4:6,10:12) = - d2Ede32 - d2Ede3de4 - d2Ede2de4 - d2Ede2de3;
dJ(10:12,4:6) = dJ(4:6,10:12)';

dJ(7:9,10:12) = - d2Ede42 - d2Ede4de3 - d2Ede1de4 - d2Ede1de3;
dJ(10:12,7:9) = dJ(7:9,10:12)';

dJ = - dJ * EI;

end
