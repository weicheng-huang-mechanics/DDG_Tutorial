% DDG tutorial, 2d_curve
% Weicheng Huang, weicheng.huang@ncl.ac.uk

clear all;
close all;
clc;

% Discrete simulation
fprintf('Beam snapping simulation \n');

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

% Compress a beam to induce buckling
totalCompress = 0.0;
inputCompress = 0.3;

% Compute first edge length
x_1 = getVertex(rodParams.x, 1);
x_2 = getVertex(rodParams.x, 2);
delta = norm(x_2 - x_1);

% Open file for writing
fileID = fopen('data.txt', 'w');  

% Simulation loop
for timeStep=1:simParams.Nsteps
    
    fprintf('t=%f\n', ctime);
    
    % Step 1: Compress the beam
    if (ctime < 3.0 && totalCompress <= inputCompress)
        
        % Add perturbation
        rodParams.g = [0.0; 10.0];
        
        rodParams.x(1) = rodParams.x(1) + simParams.dt * 0.1;
        rodParams.x(3) = rodParams.x(3) + simParams.dt * 0.1;
        rodParams.x(2*rodParams.nv-3) = rodParams.x(2*rodParams.nv-3) - simParams.dt * 0.1;
        rodParams.x(2*rodParams.nv-1) = rodParams.x(2*rodParams.nv-1) - simParams.dt * 0.1;
        
        totalCompress = totalCompress + 2 * simParams.dt * 0.1;
    end
    
    % Step 2: Remove perturbation
    if (ctime > 3.0)
        rodParams.g = [0.0; 0.0];
    end
    
    % Step 3: Rotate one end
    if (ctime > 5.0)
        x_1 = getVertex(rodParams.x, 1);
        
        edgeV = zeros(2,1);
        edgeV(1) =   1.0 * cos(  (ctime-5.0) * 0.1 );
        edgeV(2) = - 1.0 * sin(  (ctime-5.0) * 0.1 );
        
        x_2 = x_1 + delta * edgeV;
        
        rodParams.x(3) = x_2(1);
        rodParams.x(4) = x_2(2);
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
    
    % Plot figure
    if (mod(timeStep, simParams.plotStep) == 0)
        plotBeam(rodParams.x);
    end
    
    % Cout data
    if (ctime > 5.0)
        for i = 1:rodParams.nv
            xCurrent = getVertex(rodParams.x, i);
            fprintf(fileID, '%.4f %.4f %.4f \n', [ctime-5.0 xCurrent']);  % Custom formatting
        end 
    end
end

fclose(fileID);
