% Solve constraint equation
%
% -2*sq_zeta(1)/(-2*k+1)^2 + sq_zeta(2)/(k+1)^2 + 2*sq_zeta(2)/(2*k+1)^2 - sgn*1 = 0
%
function r = SolveConstraintEquation(sq_zeta, sgn)

E = [-2, 1, 2];
ALMOST_ZERO = 1e-5;

ksi1 = sq_zeta(1)*E(1);
ksi2 = sq_zeta(2)*E(2);
ksi3 = sq_zeta(3)*E(3);

TOLERANCE = 1e-12;
MAX_ITERATION = 60;

r = [];
% A. sq_zeta(1) != 0 && sq_zeta(2) != 0 && sq_zeta(3) != 0
if sq_zeta(1) > ALMOST_ZERO && sq_zeta(2) > ALMOST_ZERO && sq_zeta(3) > ALMOST_ZERO
  % First, find two roots via bisection
  if sgn > 0
   % In this case, we need to find a root < -1 and a root between -1/2 and 1/2
   % 1. Finding thew root beyond -1
   b = -3/2;
   a = -5/2;
   while RationalFunction(b, sq_zeta, E, sgn) < 0
     a = b;
     b = (-1 + b) / 2;
   end
   while RationalFunction(a, sq_zeta, E, sgn) > 0
     a = a-1;
   end
   % Find the root
   [r1, accuracy1, iterations1] = Bisection1(a, b, TOLERANCE, MAX_ITERATION, sq_zeta, E, sgn);
   r = [r, r1];
  else % if sgn < 0
   % 1. Finding thew root beyond 1/2
   a = 1.1/2;
   b = 1;
   while RationalFunction(a, sq_zeta, E, sgn) > 0
     b = a;
     a = (1/2 + a) / 2;
   end
   while RationalFunction(b, sq_zeta, E, sgn) < 0
     b = b+1;
   end
   % Find the root
   [r1, accuracy1, iterations1] = Bisection1(a, b, TOLERANCE, MAX_ITERATION, sq_zeta, E, sgn);
   r = [r, r1];
     
  end
  % 2. Finding root between -1/2 and 1/2
  a = -0.9/2; b = 0.9/2;
  while RationalFunction(a, sq_zeta, E, sgn) < 0
    a = (-1/2 + a) / 2;
  end
  while RationalFunction(b, sq_zeta, E, sgn) > 0
    b = (1/2 + b) / 2;
  end
  [r2, accuracy2, iterations2] = Bisection1(a, b, TOLERANCE, MAX_ITERATION, sq_zeta, E, sgn);
  r = [r, r2];
  
  % Computing hextic coeffcients from the following:
  %
  % (k + 1)^2*(2*k + 1) = 4 k^4 + 12 k^3 + 13 k^2 + 6 k + 1
  %
  % (2*k + 1)^2*(2*k - 1) = 16*k^4 - 8*k^2 + 1
  %
  % (k + 1)^2*(2*k - 1) = 4*k^4 + 4*k^3 - 3*k^2 - 2*k + 1
  %
  % (k + 1)^2*(2*k - 1)*(2*k + 1) = 16 k^6 + 32 k^5 + 8 k^4 - 16 k^3 - 7 k^2 + 2 k + 1
  %
  
  % Adding everything up to get the coefficients
  % 
  a0 = sgn*1 - ksi1 - ksi2 - ksi3;
  a1 = sgn*2 - ksi1*6 - ksi2*0 - ksi3*(-2);
  a2 = sgn*(-7) - ksi1*13 - ksi2*(-8) - ksi3*(-3);
  a3 = sgn*(-16) - ksi1*12 - ksi2*0 - ksi3*4;
  a4 = sgn*8 - ksi1*4 - ksi2*16 - ksi3*4;
  a5 = sgn*32;
  a6 = sgn*16;
  
  % Horner twice with kappa1 and kappa2
  a = [a0; a1; a2; a3; a4; a5; a6];
  b = Horner(a, r1);
  c = Horner(b(2:end), r2);
  
  rr = roots(c(end:-1:2));
  for i = 1:length(rr)
    if abs(imag(rr(i))) < ALMOST_ZERO, r = [r, rr(i)]; end;
  end
  
% B. sq_zeta(1) = 0
elseif sq_zeta(1) < ALMOST_ZERO && sq_zeta(2) > ALMOST_ZERO && sq_zeta(3) > ALMOST_ZERO
  % In this case the equivalent polynomial is a quartic! Good news...
  a0 = sgn*1 - ksi2 - ksi3;
  a1 = sgn*6 - ksi2*4 - ksi3*2;
  a2 = sgn*13 - ksi2*4 - ksi3;
  a3 = sgn*12;
  a4 = sgn*4;
  a = [a0, a1, a2, a3, a4];
  
  rr = roots(a(end:-1:1));
  for i = 1:length(rr)
    if abs(imag(rr(i))) < ALMOST_ZERO, r = [r, rr(i)], end;
  end

% C. sq_zeta(2) = 0
elseif sq_zeta(1) > ALMOST_ZERO && sq_zeta(2) < ALMOST_ZERO && sq_zeta(3) > ALMOST_ZERO
  % Another quartic...
  a0 = sgn - ksi1 - ksi3;
  a1 = 0 - ksi1*4 - ksi3*(-4);
  a2 = sgn*(-8) - ksi1*4 - ksi3*4;
  a3 = 0;
  a4 = sgn*16;
  a = [a0, a1, a2, a3, a4];
  
  rr = roots(a(end:-1:1));
  for i = 1:length(rr)
    if abs(imag(rr(i))) < ALMOST_ZERO, r = [r, rr(i)], end;
  end

% D. sq_zeta(3) = 0
elseif sq_zeta(1) > ALMOST_ZERO && sq_zeta(2) > ALMOST_ZERO && sq_zeta(3) < ALMOST_ZERO
  a0 = sgn - ksi1 - ksi2;
  a1 = sgn*(-2) - ksi1*2 - ksi2*(-4);
  a2 = sgn*(-3) - ksi1 - ksi2*4;
  a3 = sgn*4;
  a4 = sgn*4;
  a = [a0, a1, a2, a3, a4];
  
  rr = roots(a(end:-1:1));
  for i = 1:length(rr)
    if abs(imag(rr(i))) < ALMOST_ZERO, r = [r, rr(i)], end;
  end

% E. sq_zeta(1) = 0 and sq_zeta(2) = 0
elseif sq_zeta(1) < ALMOST_ZERO && sq_zeta(2) < ALMOST_ZERO && sq_zeta(3) > ALMOST_ZERO
  % Quadratic 
  if sgn > 0
    r = [r, (-sqrt(ksi3) - 1)/2, (+sqrt(ksi3) - 1)/2];
  end
% F. sq_zeta(1) = 0 and sq_zeta(3) = 0
elseif sq_zeta(1) < ALMOST_ZERO && sq_zeta(2) > ALMOST_ZERO && sq_zeta(3) < ALMOST_ZERO
  % Quadratic 
  if sgn > 0
    r = [r, (-sqrt(ksi2) - 1), (+sqrt(ksi2) - 1)];
  end
elseif sq_zeta(1) > ALMOST_ZERO && sq_zeta(2) < ALMOST_ZERO && sq_zeta(3) < ALMOST_ZERO
  % Quadratic 
  if sgn < 0
    r = [r, (-sqrt(-ksi1) + 1)/2, (+sqrt(-ksi1) + 1)/2 ];
  end
end

end

function f = RationalFunction(kappa, sq_zeta, e, sgn)
  f = sq_zeta(1)*e(1) / (-2*kappa + 1)^2 + sq_zeta(2)*e(2) / (kappa + 1)^2 + sq_zeta(3)*e(3) / (2*kappa + 1)^2 - sgn*1; 
end
 
function df = RationalFunctionDerivative(kappa, sq_zeta, e)
  df = -4*sq_zeta(1)*e(1)/(-2*kappa + 1)^3 - 2*sq_zeta(2)*e(2) / (kappa + 1)^3 - 4*sq_zeta(3)*e(3) / (2*kappa + 1)^3; 
end

% bisection routine
function [m, accuracy, iterations] = Bisection1(a, b, tolerance, maxIterations, sq_zeta, e, sgn) 

% INPUT
% f : Function implementation (used to evaluate the function)
% a : left initial interval bound
% b : Right initial interval bound
% tolerance: Accuracy by which the root is approximayed
% maxIterations: Maximum number of iterations.

m = (a + b) / 2;
fa = RationalFunction(a, sq_zeta, e, sgn);

% evaluate the function at m
fval = RationalFunction(m, sq_zeta, e, sgn);
i = 0; % step counter
while (abs(fval) > tolerance ) && (i < maxIterations)
    if fa > 0
        if fval > 0
            a = m;
        else 
            b = m;
        end
    else
        if fval < 0
            a = m;
        else 
            b = m;
        end
    end
    m = (a + b) / 2;
    fval = RationalFunction(m, sq_zeta, e, sgn);
    i = i + 1;
end
iterations = i; % return number of steps
accuracy = abs(fval);
if (accuracy > tolerance)
    disp('Bisection did not converge below tolerance...')
end
end