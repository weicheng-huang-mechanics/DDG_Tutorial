%% DDG tutorial, 2d_curve, Case 2: Buckling of a compressive beam
% Weicheng Huang, weicheng.huang@ncl.ac.uk
% Dezhong Tong, dezhong@umich.edu
% Zhuonan Hao, znhao@g.ucla.edu

clear all;
close all;
clc;

fprintf('Buckling of a compressive beam \n');

% Input nodes
node = importdata('inputfile/node.txt');

% Input stretching element
edge = importdata('inputfile/edge.txt');

% Input bending element
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

% Compress a beam to induce buckling
totalCompress = 0.0;
inputCompress = 0.6;

% Open file for writing
fileID = fopen('data.txt', 'w');  

% Simulation loop
for timeStep=1:simParams.Nsteps
    
    fprintf('t=%f\n', ctime);
    
    % Add perturbation
    if (ctime < 1.0)
        rodParams.g = [0.0; 0.1];
    end
    
    % Compress the beam
    if (ctime > 1.0 && totalCompress <= inputCompress)
        rodParams.x(2*rodParams.nv-3) = rodParams.x(2*rodParams.nv-3) - simParams.dt * 0.1;
        rodParams.x(2*rodParams.nv-1) = rodParams.x(2*rodParams.nv-1) - simParams.dt * 0.1;
        
        totalCompress = totalCompress + simParams.dt * 0.1;
    end
    
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
    
    if (mod(timeStep, simParams.plotStep) == 0)
        plotBeam(rodParams.x);
    end
    
    % Cout data
    for i = 1:rodParams.nv
        xCurrent = getVertex(rodParams.x, i);
        fprintf(fileID, '%.4f %.4f %.4f \n', [totalCompress xCurrent']);  % Custom formatting
    end  
    
end

fclose(fileID);
