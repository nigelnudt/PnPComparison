function [R, t, it, obj_err, img_err] = objpose(P, Qp, options)
% OBJPOSE - Object pose estimation
%   OBJPOSE(P, Qp) compute the pose (exterior orientation)
%   between the 3D point set P represented in object space 
%   and its projection Qp represented in normalized image 
%   plane. It implements the algorithm described in "Fast 
%   and Globally Convergent Pose Estimation from Video 
%   Images" by Chien-Ping Lu et. al. to appear in IEEE 
%   Transaction on Pattern Analysis and Machine intelligence
%
%   INPUTS:
%     P - 3D point set arranged in a 3xn matrix
%     Qp - 2D point set arranged in a 2xn matrix
%     options - a structure specifies certain parameters in the algorithm.
% 
%      Field name       Parameter                              Default
% 
%      OPTIONS.initR    initial guess of rotation              none
%      OPTIONS.tol      Convergence tolerance:                 1e-5 
%                       abs(new_value-old_value)/old_value<tol
%      OPTIONS.epsilon  lower bound of the objective function  1e-8
%      OPTIONS.method   'SVD' use SVD for solving rotation     'QTN'
%                       'QTN' use quaternion for solving 
%                       rotation 
%   OUTPUTS:
%     R - estimated rotation matrix
%     t - estimated translation vector
%     it - number of the iterations taken
%     obj_err - object-space error associated with the estimate
%     img_err - image-space error associated with the estimate
%

TOL = 1E-5;
EPSILON = 1E-8;
METHOD = 'SVD';

if nargin >= 3
  if isfield(options, 'tol')
    TOL = options.tol;
  end
  
  if isfield(options, 'epsilon')
    EPSILON = options.epsilon;
  end
  
  if isfield(options, 'method')
    METHOD = options.method;
  end
end
  
n = size(P,2);

% move the origin to the center of P
pbar = sum(P,2)/n;
for i = 1:n
  P(:,i) = P(:,i)-pbar;
end

Q(1:3,n) = 0;
for i = 1 : n
  Q(:,i) = [Qp(:,i);1];
end

% compute projection matrices
F(1:3,1:3,1:n) = 0;
V(1:3) = 0;
for i = 1:n
  V = Q(:,i)/Q(3,i);
  F(:,:,i) = (V*V.')/(V.'*V);
end

% compute the matrix factor required to compute t
tFactor = inv(eye(3)-sum(F,3)/n)/n;

it = 0;
if isfield(options, 'initR')  % initial guess of rotation is given
  Ri = options.initR;
  sum(1:3,1) = 0;
  for i = 1:n
    sum = sum + (F(:,:,i)-eye(3))*Ri*P(:,i);
  end
  ti = tFactor*sum;
    
  % calculate error
  Qi = xform(P, Ri, ti);
  old_err = 0;
  vec(1:3,1) = 0;
  for i = 1 : n
    vec = (eye(3)-F(:,:,i))*Qi(:,i);
    old_err = old_err + dot(vec,vec);
  end
  
else % no initial guess; use weak-perspective approximation
  % compute initial pose estimate
  [Ri, ti, Qi, old_err] = abskernel(P, Q, F, tFactor, METHOD);
  it = 1;
end

% compute next pose estimate
[Ri, ti, Qi, new_err] = abskernel(P, Qi, F, tFactor, METHOD);
it = it + 1; 

while (abs((old_err-new_err)/old_err) > TOL) & (new_err > EPSILON) 
  
  old_err = new_err;
  
  % compute the optimal estimate of R
  [Ri, ti, Qi, new_err] = abskernel(P, Qi, F, tFactor, METHOD);
  it = it + 1;

end

R = Ri;
t = ti;
obj_err = sqrt(new_err/n);

if (nargout >= 5) % calculate image-space error
  Qproj = xformproj(P, Ri, ti);
  img_err = 0;
  vec(1:3,1) = 0;
  for i = 1:n
    vec = Qproj(i)-Qp(i);
    img_err = img_err + dot(vec,vec);
  end
end
img_err = sqrt(img_err/n);

% correct possible reflection w.r.t the projection center
if t(3) < 0
  R = -R;
  t = -t;
end

% get back to original refernce frame
t = t - Ri*pbar;

% end of OBJPOSE

function [R, t, Qout, err2] = abskernel(P, Q, F, G, method)
% ABSKERNEL -  Absolute orientation kernel
%   ABSKERNEL is the function for solving the
%   intermediate absolute orientation problems
%   in the inner loop of the OI pose estimation
%   algorithm
%
%   INPUTS:
%     P - the reference point set arranged as a 3xn matrix
%     Q - the point set obtained by transforming P with
%         some pose estimate (typically the last estimate)
%     F - the array of projection matrices arranged as
%         a 3x3xn array
%     G - a matrix precomputed for calculating t
%     method - 'SVD'  -> use SVD solution for rotation
%              'QTN' -> use quaterion solution for rotation
%
%
%   OUTPUTS:
%     R - estimated rotation matrix
%     t - estimated translation vector
%     Qout - the point set obtained by transforming P with
%         newest pose estimate
%     err2 - sum of squared object-space error associated 
%         with the estimate

n = size(P,2);

for i = 1:n
  Q(:,i) = F(:,:,i)*Q(:,i);
end

% compute P' and Q'
pbar = sum(P,2)/n;
qbar = sum(Q,2)/n;
for i = 1:n
  P(:,i) = P(:,i)-pbar;
  Q(:,i) = Q(:,i)-qbar;
end

if method == 'SVD' % use SVD solution
  % compute M matrix
  M(1:3,1:3) = 0;
  for i = 1:n
    M = M+P(:,i)*Q(:,i).';
  end
  
  % calculate SVD of M
  [U,S,V] = svd(M);
  
  % compute rotation matrix R
  R = V*(U.');
elseif method == 'QTN' % use quaternion solution
  % compute M matrix
  A(1:4,1:4) = 0;
  for i = 1:n
    A = A + qmatQ([1;Q(:,i)]).'*qmatW([1;P(:,i)]);
  end
  
  % Find the largest eigenvalue of A 
  eigs_options.disp = 0;
  [V,D] = eigs(A, eye(size(A)), 1, 'LM', eigs_options);
  
  % compute rotation matrix R from the quaternion that
  % corresponds to the largest egienvalue of A
  R = quat2mat(V);
end

sum(1:3,1) = 0;
for i = 1:n
  sum = sum + F(:,:,i)*R*P(:,i);
end
t = G*sum;

Qout = xform(P, R, t);

% calculate error
err2 = 0;
vec(1:3,1) = 0;
for i = 1 : n
  vec = (eye(3)-F(:,:,i))*Qout(:,i);
  err2 = err2 + dot(vec,vec);
end

% end of ABSKERNEL



