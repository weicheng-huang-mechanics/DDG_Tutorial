function xUncons = objfun(rodParams, simParams, consParams, sElement, bElement)

% Numerical parameter
maximum_iter = simParams.maximum_iter;
tol = simParams.tol;
dt = simParams.dt;

% Free DOF
unconsInd = consParams.unconsInd;

% Mass matrix
mUncons = rodParams.m(unconsInd);
mMat = diag(mUncons);

% Damping 
viscosity = rodParams.viscosity;

% DOF vector
xUncons = rodParams.x(unconsInd);
uUncons = rodParams.u(unconsInd);

% previous position
x0Uncons = rodParams.x0(consParams.unconsInd);

% number of iterations
iter = 0;

% norm of function value 
normf = tol*10;

% Newton's method
while (normf > tol)
    
    rodParams.x(unconsInd) = xUncons;
    
    % Update stretching element 
    for i = 1:rodParams.ne
        sElement(i).nodePos_1 = getVertex(rodParams.x, sElement(i).nodeIndex(1));
        sElement(i).nodePos_2 = getVertex(rodParams.x, sElement(i).nodeIndex(2));
    end
    
    % Update bending element 
    for i =1:rodParams.nb
        bElement(i).nodePos_1 = getVertex(rodParams.x, bElement(i).nodeIndex(1));
        bElement(i).nodePos_2 = getVertex(rodParams.x, bElement(i).nodeIndex(2));
        bElement(i).nodePos_3 = getVertex(rodParams.x, bElement(i).nodeIndex(3));
        
        if (bElement(i).directSign_1 > 0)
            bElement(i).theta_1 =    getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
        else
            bElement(i).theta_1 =  - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
        end
        
        if (bElement(i).directSign_2 > 0)
            bElement(i).theta_2 =    getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
        else
            bElement(i).theta_2 =  - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
        end
        
        % update frame
        bElement(i).t_1  = (bElement(i).nodePos_2 - bElement(i).nodePos_1) / norm(bElement(i).nodePos_2 - bElement(i).nodePos_1);
        bElement(i).d_11 = parallel_transport(bElement(i).d_11_old, bElement(i).t_1_old, bElement(i).t_1);
        bElement(i).d_11 = bElement(i).d_11 / norm(bElement(i).d_11);
        bElement(i).d_12 = cross(bElement(i).t_1, bElement(i).d_11);
        bElement(i).d_12 = bElement(i).d_12 / norm(bElement(i).d_12);
        
        bElement(i).t_2  = (bElement(i).nodePos_3 - bElement(i).nodePos_2) / norm(bElement(i).nodePos_3 - bElement(i).nodePos_2);
        bElement(i).d_21 = parallel_transport(bElement(i).d_21_old, bElement(i).t_2_old, bElement(i).t_2);
        bElement(i).d_21 = bElement(i).d_21 / norm(bElement(i).d_21);
        bElement(i).d_22 = cross(bElement(i).t_2, bElement(i).d_21);
        bElement(i).d_22 = bElement(i).d_22 / norm(bElement(i).d_22);
        
        % compute reference twist
        u_temp = parallel_transport(bElement(i).d_11, bElement(i).t_1, bElement(i).t_2);
        u_temp = rotateAxisAngle(u_temp, bElement(i).t_2, bElement(i).refTwist_old);
        deltaAngle = signedAngle(u_temp, bElement(i).d_21, bElement(i).t_2);
        bElement(i).refTwist = bElement(i).refTwist_old + deltaAngle;
        
        % build material frame
        cs = cos( bElement(i).theta_1 );
        ss = sin( bElement(i).theta_1 );
        bElement(i).m_11 =   cs * bElement(i).d_11 + ss * bElement(i).d_12;
        bElement(i).m_12 = - ss * bElement(i).d_11 + cs * bElement(i).d_12;
    
        cs = cos( bElement(i).theta_2 );
        ss = sin( bElement(i).theta_2 );
        bElement(i).m_21 =   cs * bElement(i).d_21 + ss * bElement(i).d_22;
        bElement(i).m_22 = - ss * bElement(i).d_21 + cs * bElement(i).d_22;
    end
    
    % Get forces
    [Fg, Jg] = getFg(rodParams);
    [Fs, Js] = getFs(rodParams, sElement);
    [Fb, Jb] = getFb(rodParams, bElement);
    [Ft, Jt] = getFt(rodParams, bElement);
    
    Forces = (Fg + Fs + Fb + Ft);
    Forces = Forces(unconsInd);
    
    % Equation of motion
    if (rodParams.ifStatic == 1)
        f = - Forces;
    else
        f = mUncons .* ( (xUncons - x0Uncons)/dt^2 - uUncons/dt ) + viscosity * mUncons .* (xUncons - x0Uncons)/dt - Forces;
    end
    
    % Manipulate the Jacobians
    Jforces = Jg + Js + Jb + Jt;
    Jforces = Jforces(unconsInd, unconsInd);
    
    if (rodParams.ifStatic == 1)
        J = - Jforces;
    else
        J = mMat/dt^2 + viscosity * mMat/dt - Jforces;
    end
    
    % Newton's update
    xUncons = xUncons -  J\f;
    rodParams.x(unconsInd) = xUncons;
    
    % Get the norm error
    normfNew = norm(f);
    
    % Update iteration number
    iter = iter+1;
    fprintf('Iter=%d, error=%e\n', iter-1, normfNew);
    normf = normfNew;
    
    if (iter>maximum_iter)
        return
    end
end

