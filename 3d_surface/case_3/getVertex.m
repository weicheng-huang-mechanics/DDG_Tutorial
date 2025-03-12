function localNode = getVertex(x, i)

localNode = zeros(3,1);

localNode(1) = x(3*(i-1)+1);
localNode(2) = x(3*(i-1)+2);
localNode(3) = x(3*(i-1)+3);

end
