function xUncons = objfun(plateParams, simParams, consParams, sElement, bElement)
% This function solves for the equilibrium/dynamic state of a simulated system
%        using Newton's method on unconstrained degrees of freedom.
%
%   Input:
%       plateParams - struct containing plate material and state parameters
%       simParams  - struct with numerical parameters
%       consParams - struct with constraints (unconstrained DOF indices)
%       sElement   - stretching element definitions
%       bElement   - bending element definitions
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
mUncons = plateParams.m(unconsInd);
mMat = diag(mUncons);

% Damping 
viscosity = plateParams.viscosity;

% DOF vector
xUncons = plateParams.x(unconsInd);
uUncons = plateParams.u(unconsInd);

% previous position
x0Uncons = plateParams.x0(consParams.unconsInd);

% number of iterations
iter = 0;

% norm of function value 
normf = tol*10;

% Newton's method
while (normf > tol)
    
    plateParams.x(unconsInd) = xUncons;
    
    % Get forces
    [Fg, Jg] = getFg(plateParams);
    [Fs, Js] = getFs(plateParams, sElement);
    [Fb, Jb] = getFb(plateParams, bElement);
    
    Forces = (Fg + Fs + Fb);
    Forces = Forces(unconsInd);
    
    % Equation of motion
    
    if (plateParams.ifStatic == 1)
        f = - Forces;
    else
        f = mUncons .* ( (xUncons - x0Uncons)/dt^2 - uUncons/dt ) + viscosity * mUncons .* (xUncons - x0Uncons)/dt - Forces;
    end
    
    % Manipulate the Jacobians
    Jforces = Jg + Js + Jb ;
    Jforces = Jforces(unconsInd, unconsInd);
   
    if (plateParams.ifStatic == 1)
        J= - Jforces;
    else
        J= mMat/dt^2 + viscosity * mMat/dt - Jforces;
    end
    
    % Newton's update
    xUncons = xUncons -  J\f;
    plateParams.x(unconsInd) = xUncons;
    
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

