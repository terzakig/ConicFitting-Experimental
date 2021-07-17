% Solve the generalized eigenvector problem Q*x = λ*C*x
% where C = [ 0 0 -2 0 0 0; 0 1 0 0 0 0; -2 0 0 0 0 0; zeros(3, 6)]
function [Ug, Eg] = SolveEllipseGeneralizedEigenvalue(Q)
  Ug = zeros(6);
  Eg = zeros(6);
  C = [ 0  0  -2  0  0  0;...
        0  1   0  0  0  0;...
       -2  0   0  0  0  0;...
       zeros(3, 6)];
  
  % The eigen decomposition of C (nullspace dimension is 3)
  U = [-sqrt(2)/2    0    -sqrt(2)/2    0     0    0;...
           0         1         0        0     0    0;...
       -sqrt(2)/2    0     sqrt(2)/2    0     0    0;...
           0         0         0        1     0    0;...
           0         0         0        0     1    0;...
           0         0         0        0     0    1];
  E = [-2  0  0  0  0  0;...
        0  1  0  0  0  0;...
        0  0  2  0  0  0;...
        zeros(3, 6)];
        
% Convention: C = U*E*U'

N = U(:, 3:6);

% Find a 3D set of unit vectirs that don't belong to the null space of Q
%dim = 3;
%basis_found = false;
%while dim > 0 && ~basis_found
%  for i = 1:6
%    if (abs(Q*U(:, i)) < 1e-6), continue; end;
%    for j = setdiff(1:6, 1:i)
%      if (abs(Q*U(:, j)) < 1e-6), continue; end;
%        for k = 1:setdiff(1:6, union(1:i, 1:j))
%          if (abs(Q*U(:, k)) < 1e-6), continue; end;
%          H = [U(:, i), U(:, j), U(:, k)];
%          basis_found = true;
%          break;
%        end
%        if basis_found, break; end;
%    end
%    if basis_found, break; end;
%  end
%  if ~basis_found, dim = dim-1; end;      
%end
%if ~basis_found, disp('skata! no initial basis...'); continue; end;
H = U(:, 1:3);
if rank(H'*Q*H) < 3 , disp('skata! rank(H^T*Q*H) < 3'); continue; end;
W = inv(H'*C*H)*(H'*Q*H);
% Decompose W
[Uw, Ew] = eig(W);

% Now constract the remaining eigenvectors
end