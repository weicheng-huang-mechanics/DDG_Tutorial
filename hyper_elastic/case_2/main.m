% DDG tutorial, axisymmetric shell
% Weicheng Huang, weicheng.huang@ncl.ac.uk

clear all;
close all;
clc;

% Discrete Plate simulation
fprintf('Axis symmetric shell \n');

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

%
plotData = zeros(2,2);
temp = 1;

% Simulation loop
for timeStep=1:simParams.Nsteps
    
    fprintf('t=%f\n', ctime);
    
    rodParams.pressure = rodParams.pressure - 10000 * rodParams.dt;
    
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
    
    % Plot figure
    if (mod(timeStep-1, simParams.plotStep) == 0)
        plotBeam(rodParams.x);  
    end  
    
    % Cout data
    for i = 1:rodParams.nv
        xCurrent = getVertex(rodParams.x, i);
        fprintf(fileID, '%.4f %.4f %.4f \n', [ctime xCurrent']);  % Custom formatting
    end  
    
    plotData(temp,1) = rodParams.pressure;
    plotData(temp,2) = rodParams.x(1);
    temp = temp + 1;
end

fclose(fileID);


figure(2)
plot(-plotData(:,1), plotData(:,2), '-');
