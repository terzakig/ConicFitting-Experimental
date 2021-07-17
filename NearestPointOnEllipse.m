% Nearest point on ellipse to a given point
function [e, candidates] = NearestPointOnEllipse(p, u, v, c)
% p: 2D point
% u, v: Scaled axes of the ellipse), i.e. u = A*normalized(d) and v = B*normalized(v)
% c: The center of thge ellipse

% computing the solution for the ellipse stereographic parameter t
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
candidates = zeros(2, 4);
for i = 1:4
   alpha = (1-t(i)^2)/(1+t(i)^2); beta = 2*t(i)/(1+t(i)^2);
   
   m = alpha*u + beta*v + c;
   candidates(:, i) = m;
   dist = norm(m-p);
   if (dist < mindist)
     mindist = dist
     e = m
   end;
     
end 


end