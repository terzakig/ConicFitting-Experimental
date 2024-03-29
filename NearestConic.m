% Solve the nearest conic problem.
% 
% This can be stated with a Lagrangian as follows:
%
% L = (N*w + H*y)^t*(N*w + H*y) + kappa*(sign*y'*E*y + 1)
%
% where:
%
%           sign = +/-1 (ellipse/hyperbola)
%
% The null space of the constraint matrix C,
%
%           N = [ 0 0 0;...
%                 0 0 0;...
%                 0 0 0;...
%                 1 0 0;...
%                 0 1 0;...
%                 0 0 1];
%
% and the non-vanishing eigenvectors of C
%
%           H = [ -sqrt(2)/2   0   -sqrt(2)/2;...
%                      0       1        0;...
%                 -sqrt(2)/2   0    sqrt(2)/2;...
%                      0       0        0;...
%                      0       0        0;...
%                      0       0        0];
%
% The non-vanishing eigenvalues of C in a diagonal 3xs3 matrix:
%
%           E = [-2  0  0;...
%                 0  1  0;...
%                 0  0  2];


function [kappas, X] = NearestConic(v, sgn)
  N = [ 0 0 0;...
        0 0 0;...
       0 0 0;...
       1 0 0;...
       0 1 0;...
       0 0 1];

  H = [ -sqrt(2)/2   0   -sqrt(2)/2;...
            0       1        0;...
        -sqrt(2)/2   0    sqrt(2)/2;...
            0       0        0;...
            0       0        0;...
            0       0        0];

  E = [-2  0  0;...
        0  1  0;...
        0  0  2];
  
  C = [0 0 -2 0 0 0; 0 1 0 0 0 0; -2 0 0  0 0 0; zeros(3, 6)];
  
  ALMOST_ZERO = 1e-6;
  
  % Compute |Delta| = |b^2-4*a*c|
  absDelta = abs(v(2)^2 - 4*v(1)*v(3));
  % Remove scale if it applies
  if absDelta > ALMOST_ZERO, v = v / sqrt(absDelta), end;
  
  % The 3x3 vector C*H*v
  zeta = H'*v; 
  sq_zeta = zeta.^2;
  sq_zeta_norm = sum(sq_zeta);
  
  kappas = [];
  X = [];
  
  % 1. First, check if we can have solutions for k = 1/2, k = -1/2 and k = -1
  %    In any such case, dim(null(κ*C+eye())) = 1 and the only way to have a solution is
  %    if v = λ*n, where n is the null vector of κ*C+eye().
  
  % CASE #1: k = -1 and null(κ*C+eye()) = [0 -1 0 0 0 0]'
  if abs([0 -1 0 0 0 0]*v / norm(v))>0.9999
    kappas = [kappas, -1];
    X = [X, [0; -1; 0; 0 ; 0 ;0]];
  end
  
  % CASE #2: k = -1/2 and null(κ*C+eye()) = [-sqrt(2)/2 0 sqrt(2)/2 0 0 0]'
  if abs([-sqrt(2)/2 0 sqrt(2)/2 0 0 0]*v / norm(v))>0.9999
    kappas = [kappas, -1/2];
    X = [X, [-sqrt(2)/2; 0; sqrt(2)/2; 0 ; 0 ;0]];
  end
  
  % CASE #3: k = 1/2 and null(κ*C+eye()) = [sqrt(2)/2 0 sqrt(2)/2 0 0 0]'
  if abs([sqrt(2)/2 0 sqrt(2)/2 0 0 0]*v / norm(v))>0.9999
    kappas = [kappas, 1/2];
    X = [X, [sqrt(2)/2; 0; sqrt(2)/2; 0 ; 0 ;0]];
  end
  
  % Now solving for every other possible REAL value of κ
  non_trivial_kappas = SolveConstraintEquation(sq_zeta, sgn);
  for i = 1:length(non_trivial_kappas)
    k = non_trivial_kappas(i);
    kappas = [kappas, k];
    W = k*C + eye(6);
    X = [X, inv(W)*v];
  end
 
end

