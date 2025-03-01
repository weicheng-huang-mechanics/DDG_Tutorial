% DDG tutorial, multiple DOF system 
% Weicheng Huang, weicheng.huang@ncl.ac.uk

clear all;
close all;
clc;

fprintf('Multi DOF simulation \n');

% Define struct
systemParams = struct();

% Size of system, as the first node is fixed, the total DOF = N-1
systemParams.N = 4;

% Mass
systemParams.m = [0;1;1;1];
systemParams.massMatrix = diag(systemParams.m);

% viscosity
systemParams.c = [0.1;0.1;0.1]; 

% Spring stiffess
systemParams.k = [1;2;3]; 

% Initial spring length 
systemParams.l0 = [1;1;1];

% Driving amplitude
systemParams.F0 = [0.0;1.0;2.0;3.0];

% Driving frequency
systemParams.omega = [0.0;1.0;2.0;3.0]; % driving frequency

% Define struct
simParams = struct();

% Total time
simParams.totalTime = 10.0;

% Time step size
simParams.dt = 1e-2; 

% Tolerance
simParams.eps = 1e-6;

% Total step
simParams.Nsteps = round(simParams.totalTime / simParams.dt);

% Initialize 
t = linspace(0, simParams.totalTime, simParams.Nsteps);
x = zeros(systemParams.N, simParams.Nsteps); 
u = zeros(systemParams.N, simParams.Nsteps); 

% Initial conditions
x0 = [0; 1; 2; 3];
x(1, 1) = 0;
x(2, 1) = 1;
x(3, 1) = 2;
x(4, 1) = 3;
u(1, 1) = 0;
u(2, 1) = 0;
u(3, 1) = 0;
u(4, 1) = 0;

% Simulation loop
for ii=1:simParams.Nsteps-1
    
    % update DOF
    x_old = x(:,ii);
    u_old = u(:,ii);     
    
    % update time
    t_new = t(ii+1);  
    
    % Initial guess
    x_new = x_old + u_old * simParams.dt; 
    
    % Initialize error to a large value
    err = simParams.eps * 100; 
    
    totalStep = 1;
    
    % Newton's solver
    while err > simParams.eps
                
        [dF, dJ] = getF(x_new, x_old, u_old, t_new, systemParams, simParams);
        
        dF = dF(2:systemParams.N);
        dJ = dJ(2:systemParams.N,2:systemParams.N);
        
        deltaX = dJ \ dF;
        
        x_new(2:systemParams.N) = x_new(2:systemParams.N) - deltaX;
        
        err = norm(dF);
        
        totalStep = totalStep + 1;
        
        if (totalStep > 10)
            break;
        end
    end
    
    % update velocity
    u_new = (x_new - x_old) / simParams.dt;
    
    % store solution
    x(:, ii+1) = x_new; 
    u(:, ii+1) = u_new;
end

% plot result
h1 = figure(1);
plot(t, x(1,:),'k-');
hold on
plot(t, x(2,:), 'b-');
plot(t, x(3,:), 'g-');
plot(t, x(4,:), 'r-');

xlabel('t'); 
ylabel('x');
