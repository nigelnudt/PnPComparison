function R = rot2mat(L, theta)
% ROT2MAT - Converts a rotation represented by
%   rotation axis and rotation angle into a
%   3x3 orthogonal matrix
%
% INPUTS:
%   L - rotation axis, a 3x1 matrix
%   theta - rotation angle

N = L/norm(L,2);
A = [0 -N(3) N(2); N(3) 0 -N(1); -N(2) N(1) 0];
ctheta = cos(theta*pi/180);
stheta = sin(theta*pi/180);

R = ctheta*eye(3)+(1-ctheta)*N*N.'+stheta*A;
