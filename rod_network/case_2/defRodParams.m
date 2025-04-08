function rodParams = defRodParams(node, edge, bend, simParams)
% This function defines a rod struct used for the DDG-based simulation
%   Input:
%       node - nodal coordinates of the rod (nv x 2)
%       edge - edge connectivity list (ne x 2)
%       bend - bending element list (nb x 2)
%       simParams - numerical parameters of the rod
%   Output:
%       rodParams - the defined rod struct contains the physical and
%       numerical parameters of the simulated system

% Define struct
rodParams = struct();

% Young's modulus
rodParams.Y = 1e8;

% Poisson ratio
rodParams.nu = 0.5;

% Shear modulus
rodParams.G = rodParams.Y/(2.0*(1.0+rodParams.nu));

% Density
rodParams.rho = 1e3;

% Cross-sectional radius of rod
rodParams.r0 = 1e-2;

% Viscosity
rodParams.viscosity = 1.0;

% Time step size
rodParams.dt = simParams.dt;

% Material properties
rodParams.EI = rodParams.Y * pi * rodParams.r0^4 / 4;
rodParams.EA = rodParams.Y * pi * rodParams.r0^2;
rodParams.GJ = rodParams.G * pi * rodParams.r0^4/2;

% total nodal number
[rodParams.nv, ~] = size(node);

% total edge number
[rodParams.ne, ~] = size(edge);

% total bend number
[rodParams.nb, ~] = size(bend);

% total DOF size
rodParams.ndof = 3 * rodParams.nv + rodParams.ne;

% DOF vector;
rodParams.x0 = zeros(rodParams.ndof, 1);
for c=1:rodParams.nv
    rodParams.x0( 3 * (c-1) + 1) = node(c,1);
    rodParams.x0( 3 * (c-1) + 2) = node(c,2);
    rodParams.x0( 3 * (c-1) + 3) = node(c,3);
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
    
    dm = rodParams.rho * deltaLength * rodParams.r0 / 2;
    
    % mass for node
    rodParams.m(3 * (index1-1) + 1) = rodParams.m(3 * (index1-1) + 1) + dm;
    rodParams.m(3 * (index1-1) + 2) = rodParams.m(3 * (index1-1) + 2) + dm;
    rodParams.m(3 * (index1-1) + 3) = rodParams.m(3 * (index1-1) + 3) + dm;
    
    rodParams.m(3 * (index2-1) + 1) = rodParams.m(3 * (index2-1) + 1) + dm;
    rodParams.m(3 * (index2-1) + 2) = rodParams.m(3 * (index2-1) + 2) + dm;
    rodParams.m(3 * (index2-1) + 3) = rodParams.m(3 * (index2-1) + 3) + dm;
    
    % mass for rotation
    rodParams.m(3 * rodParams.nv + i) = dm * rodParams.r0^2;
end

% Gravity 
rodParams.g = [0.0; 0.0; 10.0];

% static simulation for 1, dynamic simulation for 0
rodParams.ifStatic = 0; 

end
