% The Horner scheme
function b = Horner(a, x0)

  n = length(a)-1;
  b = zeros(n+1, 1);
  b(n+1) = a(n+1);
  for i = n-1:-1:0
    b(i+1) = a(i+1) + x0*b(i+2);
  end
end