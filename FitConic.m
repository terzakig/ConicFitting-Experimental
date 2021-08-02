% Fit a conic to data points
function [X, best_error] = FitConic(points, sgn, max_iteration, tolerance, useSQP)

  % points: A 2xN or Nx2 array of points
  % sgn: -1 = ellipse, 0 = parabola, 1 = hyperbola.
  % sqp_max_iteration: Max. number of SQP iterations
  % sqp_tolerance: 
  if nargin < 5, useSQP = true; end;
  if nargin < 4
    if useSQP tolerance = 1e-5; else, tolerance = 1e-9; end;  
  end
  if nargin < 3
    if useSQP, max_iteration = 20; else, max_iteration = 30; end;
  end
  
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
    if useSQP
      if ( sgn ~= 0 && Delta*sgn < 0 ) || (Delta ~= sgn && sgn == 0)
        [e_, step] = SolveConicSQP(Q, e, sgn, max_iteration, tolerance );
      else
        e_ = e;
      end
      error = e_'*Q*e_ / (e_'*e_);
    else % Run Levenberg - Marquardt
      [U, A, B, pos] = ExtractConicParameters(e);
      [u_, v_, pos_, error] = SolveEllipseLM(A*U(:, 1), B*U(:, 2), pos, pts, tolerance, max_iteration);
      e_ = EllipseParametersToCoefs(u_, v_, pos_);
    end
    if abs(error - best_error) < ALMOST_ZERO
      X = [X, e_];
    elseif error < best_error
      X = e_;
      best_error = error;
    end
  end
  
  c = 6-n-1;
  while S(c, c) < best_error && c >0
    u = U(:, c);
    [K, X] = NearestConic(u, sgn);
    e = X(:, end);
    Delta = e(2)^2-4*e(1)*e(3);
    if useSQP
      if ( sgn ~= 0 && Delta*sgn < 0 ) || (Delta ~= sgn && sgn == 0)
        [e_, step] = SolveConicSQP(Q, e, sgn, max_iteration, tolerance );
      else
        e_ = e;
      end
      error = e_'*Q*e_ / (e_'*e_);
    else % use LM and unbiased error
      [U, A, B, pos] = ExtractConicParameters(e);
      [u_, v_, pos_, error] = SolveEllipseLM(A*U(:, 1), B*U(:, 2), pos, pts, tolerance, max_iteration);
      e_ = EllipseParametersToCoefs(u_, v_, pos_);
    end
    if abs(error - best_error) < ALMOST_ZERO
      X = [X, e_];
    elseif error < best_error
      X = e_;
      best_error = error;
    end
    c = c - 1;
  end
  
  % Now we need to translate the conic back to mu
  T = [eye(2), -mu; 0 0 1];
  for i = 1:size(X, 2)
    a = X(1, i); b = X(2, i); c = X(3, i); d = X(4, i); e = X(5, i); f = X(6, i);
    W = T'*[ a b/2 d/2; b/2 c e/2; d/2 e/2 f]*T;
    X(1, i) = W(1, 1);
    X(2, i) = W(1, 2)*2;
    X(3, i) = W(2, 2);
    X(4, i) = W(1, 3)*2;
    X(5, i) = W(2, 3)*2;
    X(6, i) = W(3, 3);
  end
end