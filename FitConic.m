% Fit a conic to data points
function [X, error] = FitConic(points, sgn)

  % points: A 2xN or Nx2 array of points
  % sgn: -1 = ellipse, 0 = parabola, 1 = hyperbola.
  
  
  ALMOST_ZERO = 1e-5;
  
  if size(points, 1) > 2, points = points', end;
  
  % Centralize the data
  mu = mean(points, 2);
  pts = points - mu;
  n = size(pts, 2);
  Q = zeros(6, 6);
  for i = 1:n
    d = [x^2; x*y; y^2; x; y; 1];
    Q = Q + d*d';
  end
  % decompose Q
  [U, S, V] = svd(Q);
  % find dimension of null space of Q
  n = 1;
  while S(9-n, 9-n) < ALMOST_ZERO
    n = n+1;
  end
  
  % Now go through the null vectors
  for i = 9:-1:10-n
    u = U(:, i);
    [K, X0] = NearestConic(u, sgn);
  end
end