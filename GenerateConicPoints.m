% Generate points on ellipse or parabola (sign = -1/+1)
function points = GenerateConicPoints(n, A, B, U, pos, std_noise, sgn)
  
  if nargin < 7, sgn = -1, end;
  if nargin < 6, std_noise = 0, end;
  if nargin < 5, pos = [0;0], end;
  if nargin < 4, U = eye(2), end;
  
  points = zeros(2, n);
  for i = 1:n
    theta = rand()*2*pi;
    x_noise = normrnd(0, std_noise);
    y_noise = normrnd(0, std_noise);
    epsilon = [x_noise; y_noise];
    if sgn < 0
      points(:, i) = pos +  U*[A*cos(theta); B*sin(theta)] + epsilon;
    elseif sgn > 0
      points(:, i) = pos +  U*[A*sec(theta); B*tan(theta)] + epsilon;
    else%if sgn == 0
      points(:, i) = pos +  U*[A*theta^2; B*2*theta] + epsilon;
    end
  end
end