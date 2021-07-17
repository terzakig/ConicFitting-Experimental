% Solve the scale invariant version of the nearest conic problem.
% 
% This can be stated with a Lagrangian as follows:
%
% L = ksi'*N*w + ksi'H*y + kappa*(y'*E*y + s) + lambda*(w'*w-1/(ksi'*ksi))
%
% where:
%
%           s = +/-1 or 0 (ellipse/hyperbola/parabola)
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


function [kappas, X] = NearestConicScaleFree(v, s)
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

    kappas = [];
    X = [];
  
  % CASE #1: Nearest parabola (s == 0)
  if s == 0
    % In this case, x = N*v / norm(N*v)
    kappas = [kappas, norm(N'*v)];
    X = [X, (N'*v) / norm(N'*v)];
  
  % CASE #2: Nearest hyperbola
  elseif s == 1
    y = -s*H'*v / norm(H'*v);
    w = -s*N'*v / norm(N'*v);
    kappas = [kappas, 1];
    x = H*y+N*w;
    X = [X, x];
  end
end