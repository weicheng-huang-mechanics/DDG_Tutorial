function plotRod(x, nv, sElement)
% This function visualizes the rod by plotting its nodal positions.
%
%   Input:
%       x - global DOF vector (ndof x 1)
%       sElement - stretch element list (ne, 2)

x1 = x(1:3:3*nv-2);
x2 = x(2:3:3*nv-1);
x3 = x(3:3:3*nv);

h1 = figure(1);

clf()
%plot3(x1,x2, x3, 'ko');

hold on

[~, ne] = size(sElement);

for i = 1:ne
    x1 = getVertex(x, sElement(i).nodeIndex(1));
    x2 = getVertex(x, sElement(i).nodeIndex(2));
    plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)], 'k-');
end

axis([-1 1 -1 1 -1 1]);
view(45,30);

drawnow

