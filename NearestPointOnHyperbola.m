% Nearest point on hyperbola to a given point
function [e, mindist, candidates] = NearestPointOnHyperbola(p, u, v, c)
% p: 2D point
% u, v: Scaled axes of the hyperbola), i.e. u = A*normalized(u) and v = B*normalized(v)
% c: The center of the hyperbola (positio corresponding to the [0; 0] in local frame)

ALMOST_ZERO = 1e-6;
% computing coefficients of the quartic that yields the possible solutions
D2 = norm(v)^2 + norm(u)^2;

X = (p - c)'*u;
Y = (p - c)'*v;

a0 = -Y;
a1 = 2*(D2-X);
a2 = 0;
a3 = 2*(D2+X);
a4 = Y;
t = roots([a4 a3 a2 a1 a0]);

% Find best solution
mindist = inf;
candidates = [];
for i = 1:4
  if abs(imag(t(i))) < ALMOST_ZERO
    alpha = (1+t(i)^2)/(1-t(i)^2); % sec
    beta = 2*t(i)/(1-t(i)^2); % tan
   
    m = alpha*u + beta*v + c;
    candidates = [candidates, m];
    dist = norm(m-p);
    if (dist < mindist)
      mindist = dist;
      e = m;
    end
  end 
end 


end