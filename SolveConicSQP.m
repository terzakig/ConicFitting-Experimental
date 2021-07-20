% SQP solution of the quadratically constrained quadratic program:
%
% f(x) = e'*Omega*e
%
% where e is a 6-vector representing a conic sectiuon
%
% s.t:
%      e'*C*e = sgn 
%      
%      C = [ 0  0 -2 0 0 0;..
%            0  1  0 0 0 0;...
%            0  0 -2 0 0 0;...
%            zeros(3, 6)]
%
%      sgn = -1 (ellipse) / 0 (parabola) / 1 (hyperbola)
%
function [e, step] = SolveConicSQP(Q, e0, sgn, maxIteration, tolerance )
%

%N = [ 0 0 0;...
%        0 0 0;...
%       0 0 0;...
%       1 0 0;...
%       0 1 0;...
%       0 0 1];

%  H = [ -sqrt(2)/2   0   -sqrt(2)/2;...
%            0       1        0;...
%        -sqrt(2)/2   0    sqrt(2)/2;...
%            0       0        0;...
%            0       0        0;...
%            0       0        0];

%  E = [-2  0  0;...
%        0  1  0;...
%        0  0  2];
  
  C = [0 0 -2 0 0 0; 0 1 0 0 0 0; -2 0 0  0 0 0; zeros(3, 6)];

  if size(e0, 1) == 1
     e0 = e0';
  end
    
  % Termination conditions
  if nargin < 5
    tolerance = 1e-5;
  end

  if nargin < 4
    maxIteration = 15;
  end

  % b. The change of the error
  delta_norm = +inf;

  % c. Timeout
  timeout = maxIteration;
  step= 1;

  % d. Loop
  while (delta_norm > tolerance) && (step < timeout)
    % Forming the linear system
    % The n-equations from the objective funtion
    W = [ Q,        2*C*e;...
           2*e'*C,     0 ];
    g = [-Q * e;...
          sgn - e'*C*e ];

     
    e_ = linsolve(w, g);
    
    delta_e = e_(1:6);
    delta_norm = norm(e_);
    
    e = e + delta_e;
    
    % increase step
    step = step + 1;
end
%E = E * sqrt(2)/ norm(e);

end
