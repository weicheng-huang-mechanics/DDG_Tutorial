function plotBeam(x)
% plotBeam Visualizes the beam (rod) by plotting its nodal positions.
%
%   Input:
%       x - Global DOF vector (2*nv x 1)

x1 = x(1:2:end);
x2 = x(2:2:end);

h1 = figure(1);
clf()
plot(x1,x2, 'k-');
hold on
axis([-0 1 -0.5 0.5]);

drawnow

