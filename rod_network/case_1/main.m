%% DDG tutorial, Net & gridshell, Case 1: Deflection of a flexible net under
% gravity
% Weicheng Huang, weicheng.huang@ncl.ac.uk
% Dezhong Tong, dezhong@umich.edu

clear all;
close all;
clc;

% Discrete simulation
fprintf('Net under gravity \n');

% input nodes
node = importdata('inputfile/node.txt');

% input stretching element
edge = importdata('inputfile/edge.txt');

% Numerical parameter
simParams = defSimParams();

% Build beam struct
rodParams = defRodParams(node, edge, simParams);

% Build stretching element
sElement = InitialStretchingElement(rodParams, edge);

% Build boundary condition
consParams = defConsParams(rodParams);

% Current time
ctime = 0.0;

% Get constrained dof and unconstrained dof
rodParams.xCons  = rodParams.x(consParams.consInd);
rodParams.xUncons = rodParams.x(consParams.unconsInd);

% Open file for writing
fileID = fopen('data.txt', 'w'); 

% Simulation loop
for timeStep=1:simParams.Nsteps
    
    fprintf('t=%f\n', ctime);
    
    % Initial guess
    rodParams.x(consParams.unconsInd) = rodParams.x(consParams.unconsInd) + rodParams.u(consParams.unconsInd) * rodParams.dt;

    % Solver
    xUncons = objfun(rodParams, simParams, consParams, sElement);
    
    % Update DOF
    rodParams.x(consParams.unconsInd) = xUncons;
    
    % Update velocity
    rodParams.u = (rodParams.x - rodParams.x0) / rodParams.dt;
    
    % Update x0
    rodParams.x0 = rodParams.x;
    
    % Update time
    ctime = ctime + simParams.dt;
    
    % Update stretching element 
    for i = 1:rodParams.ne
        sElement(i).nodePos_1 = getVertex(rodParams.x, sElement(i).nodeIndex(1));
        sElement(i).nodePos_2 = getVertex(rodParams.x, sElement(i).nodeIndex(2));
    end
    
    % Plot figure
    if (mod(timeStep, simParams.plotStep) == 0)
        plotRod(rodParams.x, rodParams.nv, sElement);
    end
    
    % Cout data
    for i = 1:rodParams.nv
        xCurrent = getVertex(rodParams.x, i);
        fprintf(fileID, '%.4f %.4f %.4f %.4f \n', [ctime xCurrent']);  % Custom formatting
    end 
end

fclose(fileID);
