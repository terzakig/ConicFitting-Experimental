% Convert Ellipse parameters to coefficients
function e_ = EllipseParametersToCoefs(u, v, pos)
  
  % u, v: The **scaled** axis direction vectoprs  (i.e., u = A*normalized(u) , v = B*normalized(v) )
  A = norm(u);
  B = norm(v);
  u1 = u(1) / A^2; u2 = u(2) / A^2;
  v1 = v(1) / B^2; v2 = v(2) / B^2;
  px = pos(1);
  py = pos(2);
  
  a = u1^2 + v1^2;
  b = 2*(u1*u2 + v1*v2);
  c = u2^2 + v2^2;
  d = -2*( (u1^2+v1^2)*px + (u1*u2+v1*v2)*py );
  e = -2*( (u2^2+v2^2)*py + (u1*u2+v1*v2)*px );
  f = (u1^2+v1^2)*px^2 + (u2^2+v2^2)*py^2 + 2*(u1*u2 + v1*v2)*px*py - 1;
  
  e_ = [a;b;c;d;e;f];
end