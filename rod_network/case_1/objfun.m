function xUncons = objfun(rodParams, simParams, consParams, sElement)
% This function solves for the equilibrium/dynamic state of a simulated system
%        using Newton's method on unconstrained degrees of freedom.
%
%   Input:
%       rodParams  - struct containing rod material and state parameters
%       simParams  - struct with numerical parameters
%       consParams - struct with constraints (unconstrained DOF indices)
%       sElement   - stretching element definitions
%
%   Output:
%       xUncons    - updated positions of unconstrained DOFs after convergence

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
    
    % Get forces
    [Fg, Jg] = getFg(rodParams);
    [Fs, Js] = getFs(rodParams, sElement);
    
    Forces = (Fg + Fs);
    Forces = Forces(unconsInd);
    
    % Equation of motion
    if (rodParams.ifStatic == 1)
        f = - Forces;
    else
        f = mUncons .* ( (xUncons - x0Uncons)/dt^2 - uUncons/dt ) + viscosity * mUncons .* (xUncons - x0Uncons)/dt - Forces;
    end
    
    % Manipulate the Jacobians
    Jforces = Jg + Js;
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

