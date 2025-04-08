function plotBeam(x)
%   Visualizes the beam by plotting its nodal positions.
%
%   Input:
%       x - global DOF vector (2*nv x 1)

x1 = x(1:2:end);
x2 = x(2:2:end);

h1 = figure(1);
clf()
plot(x1,x2, 'ko');
hold on
axis([0 0.2 -0.1 0.1]);

drawnow

