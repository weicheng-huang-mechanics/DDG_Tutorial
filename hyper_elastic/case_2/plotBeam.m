function plotBeam(x)

x1 = x(1:2:end);
x2 = x(2:2:end);

h1 = figure(1);
clf()
plot(x1,x2, 'k-');
axis([-1 9 -5 5]);

drawnow

