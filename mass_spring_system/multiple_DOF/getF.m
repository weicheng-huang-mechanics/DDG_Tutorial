function [dF, dJ] = getF(x_new, x_old, u_old, t_new, systemParams, simParams)

massMatrix = systemParams.massMatrix;
c = systemParams.c;
k = systemParams.k;
l0 = systemParams.l0;

N = systemParams.N;

F0 = systemParams.F0;
omega = systemParams.omega;

dt = simParams.dt;

dF = zeros(N,1);
dJ = zeros(N,N);

% loop over all element
for i = 1:N-1
    
    x1 = x_new(i);
    x2 = x_new(i+1);
    
    u1 = (x_new(i) - x_old(i)) / dt;
    u2 = (x_new(i+1) - x_old(i+1)) / dt;
    
    l = l0(i);
    
    % local force
    force = [- k(i) * (x2 - x1 - l) - c(i) * (u2 - u1); k(i) * (x2 - x1 - l) + c(i) * (u2 - u1) ];
    
    % local jacobian
    jac = [k(i)+c(i)/dt -k(i)-c(i)/dt; -k(i)-c(i)/dt k(i)+c(i)/dt];
    
    index(1) = i;
    index(2) = i+1;
    
    dF(index) = dF(index) + force;
    
    dJ(index,index) = dJ(index,index) + jac;
    
end

dF = massMatrix * ( (x_new - x_old) / dt - u_old ) / dt + dF - F0 .* sin(omega * t_new);

dJ = massMatrix / dt^2 - dJ;


end