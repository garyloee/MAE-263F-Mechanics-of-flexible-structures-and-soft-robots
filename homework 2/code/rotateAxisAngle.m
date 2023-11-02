function vNew = rotateAxisAngle( v, z, theta )

% This function outputs a vector "vNew" after rotating the input vector
% "v" by an angle "theta" about a vector "z".
% Example: If v = [1;0;0] (x-axis), z = [0;0;1] (z-axis), and theta=pi/2,
% then vNew = [0;1;0] (y-axis.

if (theta == 0) 
    vNew = v;
else
    c = cos(theta);
    s = sin(theta);
    vNew = c*v + s*cross(z,v) + dot(z,v) * (1.0-c) * z;
end
end
