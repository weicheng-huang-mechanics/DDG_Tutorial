function plotRod(x, nv)
% This function visualizes the rod by plotting its nodal positions.
%
%   Input:
%       x - global DOF vector (ndof x 1)

x1 = x(1:3:3*nv-2);
x2 = x(2:3:3*nv-1);
x3 = x(3:3:3*nv);

h1 = figure(1);

clf()
plot3(x1,x2, x3, 'ko-');
hold on
axis([-1 1 -1 1 -1.5 0.5]);
view(45,30);

drawnow

