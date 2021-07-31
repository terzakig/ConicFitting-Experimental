% Nearest point on patabola to a given point
function [e, mindist, candidates] = NearestPointOnParabola(p, u, v, c)
% p: 2D point
% u, v: Scaled axes of the parabola), i.e. u = A*normalized(u) and v = B*normalized(v)
% c: The "center" of the parabola (i.e., point [0; 0] in the local coordinate frame)

ALMOST_ZERO = 1e-6;
% computing coefficients of the quartic that yields the possible solutions
A2 = norm(u)^2;
B2 = norm(v)^2;

X = (p - c)'*u;
Y = (p - c)'*v;

a0 = -X;
a1 = A2 - 2*Y;
a2 = 0
a3 = 2*B2;
t = roots([a3 a2 a1 a0]);

% Find best solution
mindist = inf;
candidates = [];
for i = 1:3
  if abs(imag(t(i))) < ALMOST_ZERO
    % The parametrization convention is:
    %
    %  m = u*t + v*t^2 (i.e., in local coordinates, [A*t; B*t^2] )
    m = u*t(i) + v*t(i)^2 + c;
    candidates = [candidates, m];
    dist = norm(m-p);
    if (dist < mindist)
      mindist = dist;
      e = m;
    end
  end 
end 


end