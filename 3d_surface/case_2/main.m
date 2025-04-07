%% DDG tutorial, 3d_rod, Case 2: Plate under gravity
% Weicheng Huang, weicheng.huang@ncl.ac.uk
% Dezhong Tong, dezhong@umich.edu

clear all;
close all;
clc;

% Discrete simulation
fprintf('Plate under gravity \n');

% input nodes
node = importdata('inputfile/node.txt');

% input stretching element
edge = importdata('inputfile/edge.txt');

% input triangular element
triangle = importdata('inputfile/triangle.txt');

% Numerical parameter
simParams = defSimParams();

% Build beam struct
plateParams = defPlateParams(node, edge, triangle, simParams);

% Build stretching element
sElement = InitialStretchingElement(plateParams, edge);

% Build bending element
bElement = InitialBendingElement(plateParams, triangle);
[~, plateParams.nb] = size(bElement);

% Build boundary condition
consParams = defConsParams(plateParams);

% Current time
ctime = 0.0;

% Get constrained dof and unconstrained dof
plateParams.xCons  = plateParams.x(consParams.consInd);
plateParams.xUncons = plateParams.x(consParams.unconsInd);


% Open file for writing
fileID = fopen('data.txt', 'w'); 

% Simulation loop
for timeStep=1:simParams.Nsteps
    
    fprintf('t=%f\n', ctime);
    
    % Initial guess
    plateParams.x(consParams.unconsInd) = plateParams.x(consParams.unconsInd) + plateParams.u(consParams.unconsInd) * plateParams.dt;

    % Solver
    xUncons = objfun(plateParams, simParams, consParams, sElement, bElement);
    
    % Update DOF
    plateParams.x(consParams.unconsInd) = xUncons;
    
    % Update velocity
    plateParams.u = (plateParams.x - plateParams.x0) / plateParams.dt;
    
    % Update x0
    plateParams.x0 = plateParams.x;
    
    % Update time
    ctime = ctime + simParams.dt;
    
    % Plot figure
    if (mod(timeStep, simParams.plotStep) == 0)
        plotPlate(plateParams.x, sElement); 
    end
    
    % Cout data
    for i = 1:plateParams.nv
        xCurrent = getVertex(plateParams.x, i);
        fprintf(fileID, '%.4f %.4f %.4f %.4f \n', [ctime xCurrent']);  % Custom formatting
    end 
    
end

fclose(fileID);
