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
    
    % Get forces
    [Fg, Jg] = getFg(rodParams);
    [Fs, Js] = getFs(rodParams, sElement);
    [Fb, Jb] = getFb(rodParams, bElement);
    [Fm, Jm] = getFm(rodParams, sElement);
    [Fd, Jd] = getFd(rodParams, sElement);
    
    Forces = (Fg + Fs + Fb + Fm + Fd);
    Forces = Forces(unconsInd);
    
    % Equation of motion
    
    if (rodParams.ifStatic == 1)
        f = - Forces;
    else
        f = mUncons .* ( (xUncons - x0Uncons)/dt^2 - uUncons/dt ) - Forces;
    end
    
    % Manipulate the Jacobians
    Jforces = Jg + Js + Jb + Jm + Jd;
    Jforces = Jforces(unconsInd, unconsInd);
    
    if (rodParams.ifStatic == 1)
        J = - Jforces;
    else
        J = mMat/dt^2 - Jforces;
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

