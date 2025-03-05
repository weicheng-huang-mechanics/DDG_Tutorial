function [force, jacob] = getF(x_new, x_old, u_old, t_new, systemParams, simParams)

m = systemParams.m;
c = systemParams.c;
k = systemParams.k;

l0 = systemParams.l0;

F0 = systemParams.F0;
omega = systemParams.omega;

dt = simParams.dt;

% Equation of motion, m \ddot{x} + c \dot{x} + k x - F sin (\omega t) = 0 
force = m * ( (x_new-x_old)/dt - u_old ) / dt + c * (x_new-x_old)/dt + k * (x_new-l0) - F0 * sin(omega * t_new);

% Jacobian 
jacob = m/dt^2 + c/dt + k;

end
