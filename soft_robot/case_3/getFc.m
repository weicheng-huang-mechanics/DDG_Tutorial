function [Fc, Jc] = getFc(rodParams)
% This function computes the bending force and jacobian of the simulated 
% system. In this case the simulated system is a 2D beam
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%
% OUTPUTS: Fc - contact forces (ndof x 1)
%          Jc - contact jacobian (ndof x ndof)

Fc = zeros(rodParams.ndof, 1);
Jc = zeros(rodParams.ndof, rodParams.ndof);

dBar = 1e-3;
stiffness = 1e3;

epsilonV = 1e-6;
mu = 0.3;

for c=1:rodParams.numLeg
    
    local_node = getVertex(rodParams.x, rodParams.legEnd(c));
    
    d = local_node(2);

    if (d <= dBar)
        
        dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;
        
        Fc(2 * (rodParams.legEnd(c) - 1) + 2) = - stiffness * dEnergydD;
        
        d2EnergydD2 = - 2 * log(d / dBar) - 2 * (d - dBar) / d - 2 * (d - dBar) / d + (d - dBar) * (d - dBar) / (d * d);

        Jc(2 * (rodParams.legEnd(c) - 1) + 2, 2 * (rodParams.legEnd(c) - 1) + 2) = - stiffness * d2EnergydD2;
    end
    
    
    local_node = getVertex(rodParams.x0, rodParams.legEnd(c));
    
    d = local_node(2);
    
    if (d <= dBar)
        dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;
        
        %local_v = (rodParams.x( 2*(rodParams.legEnd(c)-1) + 1) - rodParams.x0( 2*(rodParams.legEnd(c)-1) + 1)) / rodParams.dt;
        
        local_v = rodParams.u( 2*(rodParams.legEnd(c)-1) + 1);
        
        vTangent(1) = local_v;
        vTangent(2) = 0.0;

        % the column friction is
        % 1. when v > epsilon, f_fr = mu * lam
        % 2. when v < epsilon, f_fr = f * mu * lam
        %    f = - v^2/epsilon^2 + 2 * v / epsilon
        %    f'= - 2*v/epsilon + 2/epsilon;
        %  dFf = mu * lam * ((f' * |u_t| - f)/|u_t|^3 * u_t^2 + f1 /
        % |u_t|)
        
        if ( norm(vTangent) >= epsilonV )
            f = 1.0;
            df = 0;
            tK = vTangent / norm(vTangent);
        else
            f = - ( norm(vTangent) * norm(vTangent) ) / (epsilonV * epsilonV) + 2 * norm(vTangent) / epsilonV;
            df = - 2 * norm(vTangent)/ (epsilonV * epsilonV) + 2 / epsilonV;
            tK = vTangent / (norm(vTangent) + 1e-15);
        end
        
        friction = mu * abs(stiffness * dEnergydD) * f * tK(1);
        lam = abs(stiffness * dEnergydD);
        u_t = local_v;
        dfriction = mu * lam * ((df * abs(u_t) - f)/(abs(u_t)^3) * u_t^2 + f/abs(u_t));


        if (norm(vTangent) == 0)
            friction = 0;
            dfriction = 0;
        end


        Fc(2 * (rodParams.legEnd(c) - 1) + 1) = - friction;
        
        % Jc(2 * (rodParams.legEnd(c) - 1) + 1, 2 * (rodParams.legEnd(c) - 1) + 1) = - dfriction;

    end
end

end
