function [dF, dJ] = twistingForce(node0, node1, node2, theta_e, theta_f, refTwist, ... 
    l_k, GJ, sign_1, sign_2)
% This function computes the twisting force and jacobian of a twisting
% element.
% INPUTS: node0 - position of the first node in the bending element
%         node1 - position of the second node in the bending element
%         node2 - poisiton of the third node in the bending element
%         theta_e - rotation angle of the first edge 
%         theta_f - rotation angle of the second edge
%         refTwist - the twist between the reference frames via time
%         marches
%         l_k - voronoi length of the bending element
%         GJ - twisting stiffness
%         sign_1 - the sign for the first edge
%         sing_2 - the sign for the second edge
%
% OUTPUTS: dF - twisting forces (11 x 1)
%          dJ - twisting jacobian (11 x 11)

gradTwist = zeros(11,1);

ee = node1 - node0;
ef = node2 - node1;

norm_e = norm(ee);
norm_f = norm(ef);

norm2_e = norm_e^2;
norm2_f = norm_f^2;

te = ee / norm_e;
tf = ef / norm_f;

kb = 2.0 * cross(te, tf) / (1.0 + dot(te, tf));

gradTwist(1:3) = -0.5 / norm_e * kb;
gradTwist(9:11) = 0.5 / norm_f * kb;
gradTwist(5:7) = -(gradTwist(1:3)+gradTwist(9:11));
gradTwist(4) = -1;
gradTwist(8) = 1;

if (sign_1 < 0)
    gradTwist(4) = - gradTwist(4);
end

if (sign_2 < 0)
    gradTwist(8) = - gradTwist(8);
end


DDtwist = zeros(11, 11);

chi = 1.0 + dot(te, tf);
tilde_t = (te + tf) / chi;

D2mDe2 = -0.25 / norm2_e * ( kb' * (te + tilde_t) ...
    + (te + tilde_t)' * kb);
D2mDf2 = -0.25 / norm2_f  * ( kb' * (tf + tilde_t) ...
    + (tf + tilde_t)' * kb );
D2mDeDf = 0.5 / ( norm_e * norm_f ) * ( 2.0 / chi * crossMat( te ) - ...
    kb' * tilde_t );
D2mDfDe = D2mDeDf';

DDtwist(1:3,1:3) = D2mDe2;
DDtwist(1:3, 5:7) = -D2mDe2 + D2mDeDf;
DDtwist(5:7, 1:3) = -D2mDe2 + D2mDfDe;
DDtwist(5:7, 5:7) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
DDtwist(1:3, 9:11) = -D2mDeDf;
DDtwist(9:11, 1:3) = -D2mDfDe;
DDtwist(9:11, 5:7) = D2mDfDe - D2mDf2;
DDtwist(5:7,9:11) = D2mDeDf - D2mDf2;
DDtwist(9:11,9:11) = D2mDf2;
    
integratedTwist = theta_f - theta_e + refTwist;
dF = GJ/l_k * integratedTwist * gradTwist;

dJ = GJ/l_k * (integratedTwist * DDtwist + gradTwist*gradTwist');

end
