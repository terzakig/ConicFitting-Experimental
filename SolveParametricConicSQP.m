% SQP solution of the quadratically constrained quadratic program:
%
% minimize Sum || u*x_i + v*y_i + pos - pt_i||^2
%
% subject to: u'*v = 0 
%
% where xi, yi are the coordinates of the nearest point of the conic to pti
% and u, v are the axes of the conic (rotated [1 0]' and  [0 1]' )
%
function [u, v, c, sq_error, step] = SolveParametricConicSQP(u0, v0, pos0, points, sgn, max_iteration, tolerance )
    
  % Termination conditions
  if nargin < 7
    tolerance = 1e-5;
  end

  if nargin < 6
    max_iteration = 15;
  end

  u = u0; v = v0; pos = pos0;
  theta = [u; v; pos];
  
  % Obtain the [xi, yi]
  [conic_points, sq_error] = NearestPointsOnConic(u, v, pos, sgn, points);
  
  % b. The change of the error
  delta_norm = +inf;

  % c. Timeout
  timeout = max_iteration;
  step= 1;

  % d. Loop
  C = [zeros(2), eye(2), zeros(2)]'*[eye(2), zeros(2, 4)]; % orthogonality constraint matrix
  C = C + C'; % make it a symmetric matrix
  while (delta_norm > tolerance) && (step < timeout)
    % Forming the linear system
    [Q, w] = ComputePSDMatrix(u, v, pos, conic_points, points);
    W = [ Q,             C*theta;...
          theta'*C,     0 ];
    g = [ -w;...
          -u'*v ];
     
    delta = linsolve(W, g);
    
    delta_theta = delta(1:6);
    delta_norm = norm(delta);
    
    theta = theta + delta_theta;
    u = theta(1:2);
    v = theta(3:4);
    pos = theta(5:6);
    
    % Compute the new [xi, yi] and squareed error
    [conic_points, sq_error] = NearestPointsOnConic(u, v, pos, sgn, points);
    
    % increase step
    step = step + 1;
end

end


% Compute nearest points on the conic
function [conic_points, sq_error] = NearestPointsOnConic(u, v, c, sgn, points)
  n = size(points, 2);
  conic_points = zeros(2, n);
  sq_error = 0;
  for i = 1:n
    p = points(:, i);
    if sgn < 0
      [ep, dist] = NearestPointOnEllipse(p, u, v, c);
    elseif sgn > 0
      [ep, dist] = NearestPointOnHyperbola(p, u, v, c);
    else % if sgn == 0
      [ep, dist] = NearestPointOnParabola(p, u, v, c);
    end
    conic_points(:, i) = ep;
    sq_error = sq_error + dist^2;
  end
  sq_error = sq_error / n;
end

% Compute the PSD matrix and vector of the Lagrangian from the points on the conic, i.e.
%
% L = x'*Q*x + q'*x
%
function [Q, w] = ComputePSDMatrix(u, v, c, conic_points, data_points)
  n = size(data_points, 2);
  Q = zeros(6);
  w = zeros(6, 1);
  theta = [u; v; c];
  for i = 1:n
    xi = conic_points(1, i);
    yi = conic_points(2, i);
    pi = data_points(:, i);
    
    Ji = [xi*eye(2), yi*eye(2), eye(2)];
    Q = Q + Ji'*Ji;
    w = w - Ji'*pi;
  end
    w = w + Q*theta;
    Q = Q;
end