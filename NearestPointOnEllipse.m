% Nearest point on ellipse to a given point
function [e, candidates] = NearestPointOnEllipse(p, u, v, c)
% p: 2D point
% u, v: Scaled axes of the ellipse), i.e. u = A*normalized(u) and v = B*normalized(v)
% c: The center of thge ellipse

ALMOST_ZERO = 1e-6;

% Setting up the quartic that yields up to 4 possible solutions
D2 = norm(v)^2 - norm(u)^2;

X = (p - c)'*u;
Y = (p - c)'*v;

a0 = -Y;
a1 = 2*(D2+X);
a2 = 0;
a3 = 2*(X-D2);
a4 = Y;
t = roots([a4 a3 a2 a1 a0]);

% Find best solution
mindist = inf;
candidates = [];
for i = 1:4
   if abs(imag(t(i))) < ALMOST_ZERO
    alpha = (1-t(i)^2)/(1+t(i)^2); beta = 2*t(i)/(1+t(i)^2);
   
    m = alpha*u + beta*v + c;
    candidates = [candidates, m];
    dist = norm(m-p);
    if (dist < mindist)
      mindist = dist
      e = m
    end
   end  
end 


end