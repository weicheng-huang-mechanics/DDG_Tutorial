function plotPlate(x, sElement)
% This function visualizes the plate by plotting its nodal positions.
%
%   Input:
%       x - Global DOF vector (ndof x 1)

x1 = x(1:3:end);
x2 = x(2:3:end);
x3 = x(3:3:end);

h1 = figure(1);

clf()
plot3(x1,x2, x3, 'ko');
hold on

[~, ne] = size(sElement);

for i = 1:ne
    x1 = getVertex(x, sElement(i).nodeIndex(1));
    x2 = getVertex(x, sElement(i).nodeIndex(2));
    plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)], 'k-');
end

axis([0 1 -0.5 0.5 -0.5 0.5] * 2.0);
view(25,30);

% [nv,~] = size(x);
% nv = nv / 3;
% for i = 1:50
%     plot(x1(i), x2(i), 'ko');
%     b = num2str(i);
%     c = cellstr(b);
%     text(x1(i)  + 0.01, x2(i) + 0.01, c);
% end

drawnow

