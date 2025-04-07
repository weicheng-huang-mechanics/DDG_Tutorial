%% DDG tutorial, single DOF system 
% Weicheng Huang, weicheng.huang@ncl.ac.uk
% Dezhong Tong, dezhong@umich.edu

clear all;
close all;
clc;

%% Define a struct for physical system
systemParams = struct();

% Mass
systemParams.m = 1.0;

% Spring stiffess
systemParams.k = 10; 

% Viscosity
systemParams.c = 0.1; 

% Original length of spring
systemParams.l0 = 1.0;

% Driving amplitude
systemParams.F0 = 1.0;

% Driving frequency
systemParams.omega = 1; % driving frequency

%% Define a struct for simulation parameters
simParams = struct();

% Total time
simParams.totalTime = 10;

% Time step size
simParams.dt = 1e-2; 

% Tolerance
simParams.eps = 1e-6;

% Total step
simParams.Nsteps = round(simParams.totalTime / simParams.dt);

%% Define numerical parameters
t = linspace(0, simParams.totalTime, simParams.Nsteps);
x = zeros(simParams.Nsteps, 1); 
u = zeros(simParams.Nsteps, 1); 

%% Define initial conditions
x(1) = 1.0;
u(1) = 0.0;

%% Simulation loop
for k=1:simParams.Nsteps-1 
    
    % update DOF
    x_old = x(k);
    u_old = u(k);    
    
    % update time
    t_new = t(k+1);
    
    % Initial guess
    x_new = x_old + u_old * simParams.dt; 
    
    % Initialize error to a large value
    err = simParams.eps * 100; 
    
    % Newton's solver
    while err > simParams.eps
        
        [force, jacob] = getF(x_new, x_old, u_old, t_new, systemParams, simParams);
        deltaX = force / jacob;
        
        x_new = x_new - deltaX;
        
        err = abs(force);
    end
    
    u_new = (x_new - x_old) / simParams.dt;
    
    x(k+1) = x_new; % store solution
    u(k+1) = u_new;
end

% plot result
plot(t, x, 'k-'); 
xlabel('t')
ylabel('x')
