function A = crossMat(a)
% This function computes the skew-symmetric matrix corresponding to a 3D 
% vector. e.g., A = crossMat(a) returns a 3x3 skew-symmetric matrix A such 
% that for any 3D vector b, the cross product a × b can be computed as 
% A * b.
%   Input:
%       a - a 3x1 vector
%
%   Output:
%       A - a 3x3 skew-symmetric matrix corresponding to vector a

A = [0, -a(3), a(2); ...
    a(3), 0, -a(1); ...
    -a(2), a(1), 0];
end
