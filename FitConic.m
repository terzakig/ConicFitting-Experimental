% Fit a conic to data points
function [X, best_error] = FitConic(points, sgn, sqp_max_oteration, sqp_tolerance)

  % points: A 2xN or Nx2 array of points
  % sgn: -1 = ellipse, 0 = parabola, 1 = hyperbola.
  % sqp_max_iteration: Max. number of SQP iterations
  % sqp_tolerance: 
  if nargin < 4, sqp_tolerance = 1e-5; end;
  if nargin < 3, sqp_max_iteration = 10; end;

  
  ALMOST_ZERO = 1e-5;
  
  if size(points, 1) > 2, points = points', end;
  
  % Centralize the data
  mu = mean(points, 2);
  pts = points - mu;
  n = size(pts, 2);
  Q = zeros(6, 6);
  for i = 1:n
    x = pts(1, i);
    y = pts(2, i);
    d = [x^2; x*y; y^2; x; y; 1];
    Q = Q + d*d';
  end
  % decompose Q
  [U, S, V] = svd(Q);
  % find dimension of null space of Q
  n = 1;
  while S(6-n, 6-n) < ALMOST_ZERO
    n = n+1;
  end
  
  % Now go through the null vectors
  solutions = [];
  best_error = +inf;
  for i = 6:-1:7-n
    u = U(:, i);
    [K, X] = NearestConic(u, sgn);
    
    e = X(:, end);
    Delta = u(2)^2-4*u(1)*u(3);
    if ( sgn ~= 0 && Delta*sgn < 0 ) || (Delta ~= sgn && sgn == 0)
      [e_, step] = SolveConicSQP(Q, e, sgn, sqp_max_iteration, sqp_tolerance );
    else
      e_ = e;
    end
    error = e_'*Q*e_ / (e_'*e_);
    if abs(error - best_error) < ALMOST_ZERO
      X = [X, e_];
    elseif error < best_error
      X = e_;
      best_error = error;
    end
  end
  
  c = 6-n-1;
  while S(c, c) < best_error && c >0
    Delta = e(2)^2-4*e(1)*e(3);
    if ( sgn ~= 0 && Delta*sgn < 0 ) || (Delta ~= sgn && sgn == 0)
      [e_, step] = SolveConicSQP(Q, e, sgn, sqp_max_iteration, sqp_tolerance );
    else
      e_ = e;
    end
    
    error = e_'*Q*e_ / (e_'*e_);
    if abs(error - best_error) < ALMOST_ZERO
      X = [X, e_];
    elseif error < best_error
      X = e_;
      best_error = error;
    end
    c = c - 1;
  end
end