function plateParams = defPlateParams(node, edge, triangle, simParams)
% This function defines the a plate struct used for the DDG-based simulation
%   Input:
%       node - nodal coordinates of the rod (nv x 3)
%       edge - edge connectivity list (ne x 2)
%       triangle - mesh element list (nt x 3)
%       simParams - numerical parameters of the plate
%   Output:
%       plateParams - the defined plate struct contains the physical and
%       numerical parameters of the simulated system

% Define struct
plateParams = struct();

% Young's modulus
plateParams.Y = 1e6;

% Density
plateParams.rho = 1e3;

% Plate thickness
plateParams.r0 = 5e-2;

% Viscosity
plateParams.viscosity = 0.01;

% Time step size
plateParams.dt = simParams.dt;

% Material properties
plateParams.EI = 2 * plateParams.Y * plateParams.r0^3 / (12 * sqrt(3)); 
plateParams.EA = 2 * plateParams.Y * plateParams.r0 * sqrt(3) / 4;

% total nodal number
[plateParams.nv, ~] = size(node);

% total edge number
[plateParams.ne, ~] = size(edge);

% total bend number
[plateParams.nt, ~] = size(triangle);

% total DOF size
plateParams.ndof = 3 * plateParams.nv;

% DOF vector;
plateParams.x0 = zeros(plateParams.ndof, 1);
for c=1:plateParams.nv
    plateParams.x0( 3 * (c-1) + 1) = node(c,1);
    plateParams.x0( 3 * (c-1) + 2) = node(c,2);
    plateParams.x0( 3 * (c-1) + 3) = node(c,3);
end
plateParams.x = plateParams.x0;

% Velocity vector
plateParams.u = (plateParams.x - plateParams.x0) / plateParams.dt;

% Mass array
plateParams.m = zeros(plateParams.ndof, 1);

for i=1:plateParams.nt
    
    index1 = triangle(i,1);
    index2 = triangle(i,2);
    index3 = triangle(i,3);
    
    node_1 = node(index1,:);
    node_2 = node(index2,:);
    node_3 = node(index3,:);
    
    a = norm(node_2 - node_1);
    b = norm(node_3 - node_2);
    c = norm(node_1 - node_3);
    
    s = (a+b+c) / 2;
    
    area = sqrt( s * (s-a) * (s-b) * (s-c) );
    
    dm = plateParams.rho * area * plateParams.r0 / 3;
    
    % mass for position
    plateParams.m(3 * (index1-1) + 1) = plateParams.m(3 * (index1-1) + 1) + dm;
    plateParams.m(3 * (index1-1) + 2) = plateParams.m(3 * (index1-1) + 2) + dm;
    plateParams.m(3 * (index1-1) + 3) = plateParams.m(3 * (index1-1) + 3) + dm;
    
    plateParams.m(3 * (index2-1) + 1) = plateParams.m(3 * (index2-1) + 1) + dm;
    plateParams.m(3 * (index2-1) + 2) = plateParams.m(3 * (index2-1) + 2) + dm;
    plateParams.m(3 * (index2-1) + 3) = plateParams.m(3 * (index2-1) + 3) + dm;
    
    plateParams.m(3 * (index3-1) + 1) = plateParams.m(3 * (index3-1) + 1) + dm;
    plateParams.m(3 * (index3-1) + 2) = plateParams.m(3 * (index3-1) + 2) + dm;
    plateParams.m(3 * (index3-1) + 3) = plateParams.m(3 * (index3-1) + 3) + dm;
end

% Gravity 
g = [0.0; 0.0; 0.0];

plateParams.garr = zeros(plateParams.ndof, 1);
for c = 1:plateParams.nv
    plateParams.garr( 3 * (c-1) + 1 : 3 * (c-1) + 3) = g;
end

% External force
plateParams.load = 0.0;

end
