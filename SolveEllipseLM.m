% Fit ellipse to data using the Levenberg - Marquadrt algorithm

function [u, v, c] = SolveEllipseLM(u0, v0, c0, points, tolerance, max_iteration)
  
  if nargin < 6, max_iteration = 30; end;
  if nargin < 5, tolerance = 1e-5; end;
    
  if size(points, 1) > 2, points = points'; end;  
  
  % s is applied to enforce that v is proportional to [-u0(2); u0(1)]
  % this favours the parametrization U = [u1 -u2; u2 u1]
  s = -sign(u0(2)*v0(1)); 
  v = s*v0;
  u = u0;
  U = [u, v];
  
  [sq_error, points_on_ellipse] = SquaredError(u, s*v, c, points);
  stop = sq_error < tolerance;
  lambda = 0.001;
  step = 0;
  
  while ~stop
    [Omega, ksi] = FisherParameters(u, v, s, c, points, points_on_ellipse);
    
    improvement = false;
    while ~improvement && ~stop
      x_temp = inv(Omega + lambda*eye(4))*ksi;
      u_temp = x_temp(1:2);
      v_temp = x_temp(3:4);
      c_temp = x_temp(5:6);
      [sq_error_temp, temp_points_on_ellipse] = SquaredError(u, s*v, c, points);
      if sq_error_temp < sq_error
        u = u_temp;
        v = v_temp;
        c = c_temp;
        
        lambda = lambda / 10;
        improvement = true;
      else
        lambda = lambda * 10;
      end
      step = step + 1;
      if step > max_iteration || sq_error < tolerance
        stop = true;
      end
    end
  end
  
end

% Compute nearest poinst and squared error
function [sq_error, ellipse_points] = SquaredError(u, v, c, points)
  n = size(points, 2);
  ellipse_points = zeros(2, n);
  sq_error = 0;
  for i = 1:n
    p = points(:, i);
    [ep, dist] = NearestPointOnEllipse(p, u, v, c);
    ellipse_points(:, i) = ep;
    sq_error = sq_error + dist^2;
  end
  sq_error = sq_error / n;
end

% Return Fihser partamers for a given conic and points
function [Omega, ksi] = FisherParameters(u, v, s, c, data_points, points_on_ellipse)
  n = size(points, 2);
  
  Omega = zeros(4);
  ksi = zeros(4, 1);
  
  for i = 1:n
    % Compute Jacobian
    ellipse_point = points_on_ellipse(:, i);
    data_point = data_points(:, i);
    [Df, f] = EllipseJacobian(u, v, s, c, ellipse_point, data_point);
    Omega = Omega + Df'*Df;
    ksi = ksi + Df'*f;
  end
end

% Jacobian of distance from point to (nearest point on) ellipse
function [Df, f] = EllipseJacobian(u, v, s, c, ellipse_point, data_point)
  x = ellipse_point(1);
  y = ellipse_point(2);
  Df = [x   -s*y 1 0;...
        s*y   x  0 1];
  f = u*x + s*v*y + c- data_point;
end