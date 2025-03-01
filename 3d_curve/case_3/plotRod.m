function plotRod(rodParams, bElement)

h1 = figure(1);
clf()

% define a width
widht = 0.1;

% edge
for i = 1:rodParams.nb
    localNode = bElement(i).nodePos_2;
    m1 = (bElement(i).m_11 + bElement(i).m_21) / 2;
    
    xC(i, :) = localNode';
    xL(i, :) = localNode' + widht * m1' / 2;
    xR(i, :) = localNode' - widht * m1' / 2;
end

xC(end+1,:) = xC(1,:);
xR(end+1,:) = xR(1,:);
xL(end+1,:) = xL(1,:);

plot3(xC(:,1), xC(:,2), xC(:,3), 'k-');
hold on
plot3(xL(:,1),xL(:,2), xL(:,3), 'k-');
plot3(xR(:,1),xR(:,2), xR(:,3), 'k-');

%plot3([xL(1,1) xR(1,1)],[xL(1,2) xR(1,2)],[xL(1,3) xR(1,3)],'k-');
%plot3([xL(end,1) xR(end,1)],[xL(end,2) xR(end,2)],[xL(end,3) xR(end,3)],'k-')

axis([-1 1 -1 1 -1 1]*1.2);
view(45,30);

drawnow

