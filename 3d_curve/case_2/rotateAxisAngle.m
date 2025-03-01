function u = rotateAxisAngle(v, t, angle)

if (angle == 0)
    u = v;
else
    cs = cos(angle);
	ss = sin(angle);
	u = cs*v + ss*cross(t,v)+dot(t,v)*(1.0-cs)*t;
end
    
end
