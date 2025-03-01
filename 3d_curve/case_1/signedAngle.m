function angle = signedAngle( u, v, n )

w = cross(u,v);
angle = atan2( norm(w), dot(u,v) );

if (dot(n,w) < 0) 
    angle = -angle;
end

end