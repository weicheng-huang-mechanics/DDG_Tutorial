function rodParams = defRodParams(node, edge, bend, simParams)

% This function defines a rotational shell struct used for the DDG-based simulation
%   Input:
%       node - Nodal coordinates of the rod (nv x 2)
%       edge - Edge connectivity list (ne x 2)
%       bend - Bending element list (nb x 2)
%       simParams - numerical parameters of the shell
%   Output:
%       rodParams - the defined rod struct contains the physical and
%       numerical parameters of the simulated system

% Define struct
rodParams = struct();

% Young's modulus
rodParams.Y = 1e6;

% Poisson ratio
rodParams.nu = 0.5;

% Density
rodParams.rho = 1e3;

% Thickness of thickness
rodParams.r0 = 2e-2;

% Viscosity
rodParams.viscosity = 0.1;

% Time step size
rodParams.dt = simParams.dt;

% Material properties
rodParams.EI = rodParams.Y * rodParams.r0^3 / ( 12 * (1 - rodParams.nu * rodParams.nu) );
rodParams.EA = rodParams.Y * rodParams.r0 / ( 1 - rodParams.nu * rodParams.nu );

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
    
    dm = rodParams.rho * deltaLength * pi * (node_1(1) + node_2(1)) * rodParams.r0;
    
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

% External pressure
rodParams.pressure = 0.0;

% static simulation for 1, dynamic simulation for 0
rodParams.ifStatic = 1; 

end
