function rodParams = defRodParams(node, edge, bend, simParams)
% This function defines a beam struct used for the DDG-based simulation
%   Input:
%       node - Nodal coordinates of the rod (nv x 2)
%       edge - Edge connectivity list (ne x 2)
%       bend - Bending element list (nb x 2)
%       simParams - numerical parameters of the beam
%   Output:
%       rodParams - the defined rod struct contains the physical and
%       numerical parameters of the simulated system

% Define struct
rodParams = struct();

% Young's modulus
rodParams.Y = 1e8;

% Density
rodParams.rho = 1e1;

% Cross-sectional radius of rod
rodParams.r0 = 1e-2;

% Viscosity
rodParams.viscosity = 0.1;

% Time step size
rodParams.dt = simParams.dt;

% Material properties
rodParams.EI = rodParams.Y * pi * rodParams.r0^4 / 4;
rodParams.EA = rodParams.Y * pi * rodParams.r0^2;

rodParams.crossSection = pi * rodParams.r0^2;

% total nodal number
[rodParams.nv, ~] = size(node);

% total edge number
[rodParams.ne, ~] = size(edge);

% total bend number
[rodParams.nb, ~] = size(bend);

% total DOF size
rodParams.ndof = 2 * rodParams.nv;

% DOF vector;
rodParams.x0 = zeros(rodParams.ndof, 1);
for c=1:rodParams.nv
    rodParams.x0( 2 * (c-1) + 1) = node(c,1);
    rodParams.x0( 2 * (c-1) + 2) = node(c,2);
end
rodParams.x = rodParams.x0;

% Velocity vector
rodParams.u = (rodParams.x - rodParams.x0) / rodParams.dt;

% Mass array
rodParams.m = zeros(rodParams.ndof, 1);

for i = 1:rodParams.ne
    
    index1 = edge(i,1);
    index2 = edge(i,2);
    
    node_1 = node(index1,:);
    node_2 = node(index2,:);
    
    deltaLength = norm(node_2 - node_1);
    
    dm = rodParams.rho * deltaLength * pi * rodParams.r0^2 / 2;
    
    rodParams.m(2 * (index1-1) + 1) = rodParams.m(2 * (index1-1) + 1) + dm;
    rodParams.m(2 * (index1-1) + 2) = rodParams.m(2 * (index1-1) + 2) + dm;
    
    rodParams.m(2 * (index2-1) + 1) = rodParams.m(2 * (index2-1) + 1) + dm;
    rodParams.m(2 * (index2-1) + 2) = rodParams.m(2 * (index2-1) + 2) + dm;
end

% Gravity 
g = [0.0; 0.0];

rodParams.garr = zeros(rodParams.ndof, 1);
for c = 1:rodParams.nv
    rodParams.garr( 2 * (c-1) + 1 : 2 * (c-1) + 2) = g;
end

% static simulation for 1, dynamic simulation for 0
rodParams.ifStatic = 0; 

% external magnetic field
rodParams.Ba = [1;0];

end
