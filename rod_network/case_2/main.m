% DDG tutorial, 3D rod network
% Weicheng Huang, weicheng.huang@ncl.ac.uk

clear all;
close all;
clc;

% Discrete simulation
fprintf('Buckling of gridshell \n');

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

% read target for footpoint
targetPos = importdata('inputfile/target.txt');
[targetNum, ~] = size(targetPos);

% Open file for writing
fileID = fopen('data.txt', 'w'); 

% Simulation loop
for timeStep=1:simParams.Nsteps
    
    fprintf('t=%f\n', ctime);
    
    % Add the small perturbation 
    if (ctime < 1.0)
        rodParams.g = [0.0; 0.0; 10.0];
    end
    
    % Compress the foot print
    if (ctime > 1.0)
        
        node_s = zeros(3,1);
        node_e = zeros(3,1);
        
        for kk = 1:targetNum
            
            % for strat point
            node_s(1) = targetPos(kk, 1); 
            node_s(2) = targetPos(kk, 2); 
            node_s(3) = 0.0;
            
            node_current = getVertex(rodParams.x, consParams.fixIndex (2 * (kk-1) + 1) );
            
            if ( norm(node_current - node_s) < 1e-3 )
                node_current = node_s;
            else
                tangent = (node_s - node_current) / norm(node_s - node_current);
                node_current = node_current + tangent * simParams.dt * 0.1;
            end
            
            rodParams.x( 3 * (consParams.fixIndex (2 * (kk-1) + 1)-1) + 1 ) = node_current(1);
            rodParams.x( 3 * (consParams.fixIndex (2 * (kk-1) + 1)-1) + 2 ) = node_current(2);
            
            % for end point
            node_s(1) = targetPos(kk, 3); 
            node_s(2) = targetPos(kk, 4); 
            node_s(3) = 0.0;
            
            node_current = getVertex(rodParams.x, consParams.fixIndex (2 * (kk-1) + 2) );
            
            if ( norm(node_current - node_s) < 1e-3 )
                node_current = node_s;
            else
                tangent = (node_s - node_current) / norm(node_s - node_current);
                node_current = node_current + tangent * simParams.dt * 0.1;
            end
            
            rodParams.x( 3 * (consParams.fixIndex (2 * (kk-1) + 2)-1) + 1 ) = node_current(1);
            rodParams.x( 3 * (consParams.fixIndex (2 * (kk-1) + 2)-1) + 2 ) = node_current(2);
            
        end
        
    end
    
    % Delete the small perturbation 
    if (ctime > 5.0)
        rodParams.g = [0.0; 0.0; 0.0];
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
    
    if (mod(timeStep, simParams.plotStep) == 0)
        plotRod(rodParams.x, rodParams.nv, sElement);
    end
    
    for i = 1:rodParams.nv
        xCurrent = getVertex(rodParams.x, i);
        fprintf(fileID, '%.4f %.4f %.4f %.4f \n', [ctime xCurrent']);  % Custom formatting
    end 
    
end

fclose(fileID);
