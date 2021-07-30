% Generate points on ellipse or parabola (sign = -1/+1)
function points = GenerateConicPoints(n, A, B, U, pos, std_noise, sgn)
  
  if nargin < 7, sgn = -1, end;
  if nargin < 6, std_noise = 0, end;
  if nargin < 5, pos = [0;0], end;
  if nargin < 4, U = eye(2), end;
  
  points = zeros(2, n);
  for i = 1:n
    if sgn ~= 0, theta = rand()*2*pi, else theta = rand()*10 -5, end;
      
    x_noise = normrnd(0, std_noise);
    y_noise = normrnd(0, std_noise);
    epsilon = [x_noise; y_noise];
    if sgn < 0
      points(:, i) = pos +  U*[A*cos(theta); B*sin(theta)] + epsilon;
    elseif sgn > 0
      points(:, i) = pos +  U*[A*sec(theta); B*tan(theta)] + epsilon;
    else%if sgn == 0
      points(:, i) = pos +  U*[A*theta; B*theta^2] + epsilon;
    end
  end
end
