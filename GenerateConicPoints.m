% Generate points on ellipse or parabola (sign = -1/+1)
function points = GenerateConicPoints(n, D, U, pos, std_noise, sgn)
  
  if nargin < 6, sgn = -1, end;
  if nargin < 5, std_noise = 0, end;
  if nargin < 4, pos = [0;0], end;
  if nargin < 3, U = eye(2), end;
  
  points = zeros(2, n);
  for i = 1:n
    theta = rand()*2*pi;
    x_noise = normrnd(0, std_noise);
    y_noise = normrnd(0, std_noise);
    epsilon = [x_noise; y_noise];
    if sgn == -1
      points(:, i) = pos +  D*U*[cos(theta); sin(theta)] + epsilon;
    elseif sgn = 1
      points(:, i) = pos +  D*U*[sec(theta); tan(theta)] + epsilon;
    elseif sgn == 0
      points(:, i) = pos +  D*U*[theta^2; 2*theta] + epsilon;
    end
  end
end