function localNode = getVertex(x, i)

localNode = zeros(2,1);

localNode(1) = x(2*(i-1)+1);
localNode(2) = x(2*(i-1)+2);

end
