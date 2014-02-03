function [L, theta] = mat2rot(R)
% MAT2ROT - Finds rotation axis and rotation 
%   angle from a 3x3 orthogonal matrix
%
% OUTPUTS:
%   L - rotation axis, a 3x1 matrix
%   theta - rotation angle

ctheta = (trace(R)-1)/2;
theta = acos(ctheta)*180/pi;

[V,D] = eig(R);
L(1:3,1) = 0;
for i= 1:3
  if (D(i,i) == 1)
    L = V(:,i);
    break;
  end
end

L = L/norm(L,2);


