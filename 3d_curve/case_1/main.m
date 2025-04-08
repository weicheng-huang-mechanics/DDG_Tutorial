%% DDG tutorial, 3d_rod, Case 1: Deformation of a 3D helix under gravity
% Weicheng Huang, weicheng.huang@ncl.ac.uk
% Dezhong Tong, dezhong@umich.edu
% Zhuonan Hao, znhao@g.ucla.edu

clear all;
close all;
clc;

% Discrete rod simulation
fprintf('Deformation of a 3D helix under gravity \n');

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

% Open file for writing
fileID = fopen('data.txt', 'w'); 

% Simulation loop
for timeStep=1:simParams.Nsteps
    
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
    
    % Update stretching element 
    for i = 1:rodParams.ne
        sElement(i).nodePos_1 = getVertex(rodParams.x, sElement(i).nodeIndex(1));
        sElement(i).nodePos_2 = getVertex(rodParams.x, sElement(i).nodeIndex(2));
    end
    
    % Update bending element 
    for i =1:rodParams.nb
        bElement(i).nodePos_1 = getVertex(rodParams.x, bElement(i).nodeIndex(1));
        bElement(i).nodePos_2 = getVertex(rodParams.x, bElement(i).nodeIndex(2));
        bElement(i).nodePos_3 = getVertex(rodParams.x, bElement(i).nodeIndex(3));
        
        if (bElement(i).directSign_1 > 0)
            bElement(i).theta_1 =    getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
        else
            bElement(i).theta_1 =  - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
        end
        
        if (bElement(i).directSign_2 > 0)
            bElement(i).theta_2 =    getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
        else
            bElement(i).theta_2 =  - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
        end
        
        bElement(i).t_1_old  = bElement(i).t_1;
        bElement(i).d_11_old = bElement(i).d_11;
        bElement(i).d_12_old = bElement(i).d_12; 
    
        bElement(i).t_2_old  = bElement(i).t_2;
        bElement(i).d_21_old = bElement(i).d_21;
        bElement(i).d_22_old = bElement(i).d_22;
        
        bElement(i).refTwist_old = bElement(i).refTwist;
        
        cs = cos( bElement(i).theta_1 );
        ss = sin( bElement(i).theta_1 );
        bElement(i).m_11 =   cs * bElement(i).d_11 + ss * bElement(i).d_12;
        bElement(i).m_12 = - ss * bElement(i).d_11 + cs * bElement(i).d_12;
    
        cs = cos( bElement(i).theta_2 );
        ss = sin( bElement(i).theta_2 );
        bElement(i).m_21 =   cs * bElement(i).d_21 + ss * bElement(i).d_22;
        bElement(i).m_22 = - ss * bElement(i).d_21 + cs * bElement(i).d_22;
    end
    
    % Plot figure
    if (mod(timeStep, simParams.plotStep) == 0)
        plotRod(rodParams.x, rodParams.nv);
    end
    
    % Cout data
    for i = 1:rodParams.nv
        xCurrent = getVertex(rodParams.x, i);
        fprintf(fileID, '%.4f %.4f %.4f %.4f \n', [ctime xCurrent']);  % Custom formatting
    end  
    
end

fclose(fileID);
