function plotBeam(x)
%   Visualizes the membrane by plotting its nodal positions.
%
%   Input:
%       x - global DOF vector (2*nv x 1)

x1 = x(1:2:end);
x2 = x(2:2:end);

x1 = [x1;x1(1)];
x2 = [x2;x2(1)];

h1 = figure(1);
clf()
plot(x1,x2, 'k-');
axis([-1 9 -5 5]);

drawnow

