function Rs = randrotmat(n)
% RANDROTMAT - Random rotation matrix
%   RANDROTMAT(n) generates n uniformly distributed
%   rotation matrix

Rs(1:3,1:3,1:n) = 0;
for i = 1:n
  x = rand(1,3);
  theta = 2*pi*x(1);
  phi = 2*pi*x(2);
  z = x(3);
  
  v = [cos(phi)*sqrt(z); sin(phi)*sqrt(z); sqrt(1-z)];
  
  c = cos(theta); s = sin(theta);
  Rs(:,:,i) = (2*v*v.'-eye(3))* [c,s,0;-s, c,0;0,0,1];
end
     
      


