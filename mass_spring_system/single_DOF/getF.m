function [force, jacob] = getF(x_new, x_old, u_old, t_new, systemParams, simParams)
% This function computes the force and jacobian of the simulated system. In
% this case the simulated system is a S-DOF mass-spring-damper system
% INPUTS: x_new - nodal position at the current time step.
%         x_old - nodal position at the previous time step.
%         u_old - nodal velocity at the previous time step.
%         t_new - current time stamp
%         systemParams - the physical system struct
%         simParams - the simulation parameters struct
%
% OUTPUTS: force - the sum of the forces for the simulated system.
%          jacob - the jacobian of the simulated system.

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
