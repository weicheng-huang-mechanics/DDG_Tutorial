% DDG tutorial, 2d_curve
% Weicheng Huang, weicheng.huang@ncl.ac.uk

clear all;
close all;
clc;

% Discrete beam simulation
fprintf('Beam under magnetic field \n');

% input nodes
node = importdata('inputfile/node.txt');

% input stretching element
edge = importdata('inputfile/edge.txt');

% input bending element
bend = importdata('inputfile/bend.txt');

% Numerical parameter
simParams = defSimParams();

% Build beam struct
rodParams = defRodParams(node, edge, bend, simParams);

% Build stretching element
sElement = InitialStretchingElement(rodParams, edge);

% Build bending element
bElement = InitialBendingElement(rodParams, bend, sElement);

% Build boundary condition
consParams = defConsParams(rodParams);

% Current time
ctime = 0.0;

% Get constrained dof and unconstrained dof
rodParams.xCons  = rodParams.x(consParams.consInd);
rodParams.xUncons = rodParams.x(consParams.unconsInd);

% Open file for writing
fileID = fopen('data.txt', 'w');  

% Magnetic frequency
omega = 1.0;

% Simulation loop
for timeStep=1:simParams.Nsteps
    
    
    % Update Ba
    rodParams.Ba(1) = 3 * cos(omega * ctime); 
    rodParams.Ba(2) = 3 * sin(omega * ctime); 
    
    fprintf('t=%f\n', ctime);
    
    % Initial guess
    rodParams.x(consParams.unconsInd) = rodParams.x(consParams.unconsInd) + rodParams.u(consParams.unconsInd) * rodParams.dt;

    % Solver
    xUncons = objfun(rodParams, simParams, consParams, sElement, bElement);
    
    % Update DOF
    rodParams.x(consParams.unconsInd) = xUncons;
    
    % Update velocity
    rodParams.u = (rodParams.x - rodParams.x0) / rodParams.dt;
    
    % Update x0
    rodParams.x0 = rodParams.x;
    
    % Update time
    ctime = ctime + simParams.dt;
    
    % Plot figure
    if (mod(timeStep-1, simParams.plotStep) == 0)
        plotBeam(rodParams.x);  
    end  
    
    % Cout data
    for i = 1:rodParams.nv
        xCurrent = getVertex(rodParams.x, i);
        fprintf(fileID, '%.4f %.4f %.4f \n', [ctime xCurrent']);  % Custom formatting
    end  
end

fclose(fileID);
