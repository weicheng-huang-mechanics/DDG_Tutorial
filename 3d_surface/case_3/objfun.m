function [xUncons, externalForce] = objfun(plateParams, simParams, consParams, sElement, bElement)

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
    [Fs, Js] = getFs(plateParams, sElement);
    [Fb, Jb] = getFb(plateParams, bElement);
    
    Forces = (Fs + Fb);
    
    externalForce = Forces((417-1) * 3 + 3);
    
    Forces = Forces(unconsInd);
    
    % Equation of motion
    %f = mUncons .* ( (xUncons - x0Uncons)/dt^2 - uUncons/dt ) + viscosity * mUncons .* (xUncons - x0Uncons)/dt - Forces;
    f = - Forces;
    
    % Manipulate the Jacobians
    Jforces = Js + Jb;
    Jforces = Jforces(unconsInd, unconsInd);
    
    %J= mMat/dt^2 + viscosity * mMat/dt - Jforces;
    J =- Jforces;
    
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

