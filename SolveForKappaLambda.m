% Solve the polynomial for k give a lambda, such that
%
% norm(inv(E)*Un'*ksi) = 1
%
% where:
%
% E = diag([lambda, lambda, lambda, -2*kappa+lambda, kappa+lambda, 2*kappa+lambda])
%
% and Un = ⎡         √2     -√2 ⎤
%          ⎢0  0  0  ──  0  ────⎥
%          ⎢         2       2  ⎥
%          ⎢                    ⎥
%          ⎢0  0  0  0   1   0  ⎥
%          ⎢                    ⎥
%          ⎢         √2      √2 ⎥
%          ⎢0  0  0  ──  0   ── ⎥
%          ⎢         2       2  ⎥
%          ⎢                    ⎥
%          ⎢1  0  0  0   0   0  ⎥
%          ⎢                    ⎥
%          ⎢0  1  0  0   0   0  ⎥
%          ⎢                    ⎥
%          ⎣0  0  1  0   0   0  ⎦

function [kappas, lambda, Un] = SolveForKappaLambda(ksi)
 
  % Normalize ksi
  ksi = ksi / norm(ksi);
  Un = [0   0   0   sqrt(2)/2  0  -sqrt(2)/2; ...
        0   0   0      0       1       0; ...
        0   0   0   sqrt(2)/2  0   sqrt(2)/2;...
        1   0   0      0       0       0;...
        0   1   0      0       0       0;...
        0   0   1      0       0       0];
  % Compute zeta = Un'*ksi
  zeta = Un'*ksi;
  sq_zeta = zeta.^2;
  % Now compute lambda, such that k=lambda will be one of the solutions
  lambda = sqrt( sum(sq_zeta(1:4)) + sq_zeta(5)/4 + sq_zeta(6)/9 );
  
  % Finally, the hardest part: Find the coefficients of the equivalent quartic equation, of:
  %
  %  gamma = sq_zeta(4)/(-2*kappa+lambda)^2+sq_zeta(5)/(kappa+lambda)^2+sq_zeta(6)/(2*kappa+lambda)^2
  %
  % where gamma = 1 - sum(sq_zeta(1:3))/lambda^2
  %
  gamma = 1 - sum(sq_zeta(1:3))/lambda^2;
  
  % Store the lambda-rppt as kappa0
  
  kappa0 = lambda;
  
  % Now we need to compute the root beyond -lambda with the N-R method 
  % (using -3/2*lambda as starting point)
  step = 0;
  kappa1 = -3/2*lambda;
  f1 = +inf;
  while abs(f1) > 1e-6 && step < 20
    
    f1 = RationalFunction(kappa1, sq_zeta, lambda, gamma);
    df1 = RationalFunctionDerivative(kappa1, sq_zeta, lambda);
 
    kappa1 = kappa1 - f1 / df1;
    
    step = step + 1;
  end
  
  % Great! Now we have two roots of the equivalent sextic polynomial,
  % and can therefore use the Horner scheme to get the quartic quotient
  % and then solve for the remaining roots analytically.
  
  %
  % The overall polynomial is a sum of the following:
  %
  %                      2          2                     4      3        2  2        3    4
  %  sq_zeta(4) ⋅ (k + l) ⋅ (2⋅k + l) = sq_zeta(4) ⋅ ( 4⋅k  + 4⋅k ⋅l - 3⋅k ⋅l  - 2⋅k⋅l  + l )
  
  %                        2          2                     4      2  2    4
  %  sq_zeta(5) ⋅ (2⋅k - l) ⋅(2⋅k + l) = sq_zeta(5) ⋅ ( 16⋅k  - 8⋅k ⋅l  + l )
  %
  %                      2           2                    4       3         2  2        3    4
  %  sq_zeta(6) ⋅ (k + l) ⋅ (2⋅k - l) = sq_zeta(6) ⋅ ( 4⋅k  + 12⋅k ⋅l + 13⋅k ⋅l  + 6⋅k⋅l  + l )
  %
  %                  2           2           2                  6         5     2 4       3 3      4 2      5      6
  %  -gamma ⋅ (k + l) ⋅ (2⋅k - l) ⋅ (2⋅k + l)  = -gamma ⋅ ( 16⋅k  + 32⋅l⋅k + 8⋅l⋅k  - 16⋅l⋅k  - 7⋅l⋅k  + 2⋅l⋅k  + l )
  % 
  a0 = sq_zeta(4)*lambda^4       + sq_zeta(5)*lambda^4 +     sq_zeta(6)*lambda^4     + (-gamma)*lambda^6;
  a1 = sq_zeta(4)*(-2*lambda^3)  +                           sq_zeta(6)*(6*lambda^3) + (-gamma)*(2*lambda^5);
  a2 = sq_zeta(4)*(-3*lambda^2)  + sq_zeta(5)*(-8*lambda^2) + sq_zeta(6)*(13*lambda^2) + (-gamma)*(-7*lambda^4);
  a3 = sq_zeta(4)*(4*lambda)     +                           sq_zeta(6)*(12*lambda)    + (-gamma)*(-16*lambda^3);
  a4 = sq_zeta(4)*4              + sq_zeta(5)*16            + sq_zeta(6)*4            +(-gamma)*(8*lambda^2);
  a5 =                                                                                 (-gamma)*(32*lambda);
  a6 =                                                                                 (-gamma)*16;
  
  a = [a0; a1; a2; a3; a4; a5; a6];
  
  % We now have to use the Horner scheme to get the quartic quotient polynomial
  % by dividing with (x-kappa0)*(x-kappa1)
  b = Horner(a, kappa0);
  q = Horner(b(2:end), kappa1);
  
  % So, now, q is the deflated a into a quartic, which we need to solve analytically
  
  
 end
 
 
 function f = RationalFunction(kappa, sq_zeta, lambda, gamma)
  f = sq_zeta(4) / (-2*kappa + lambda)^2 + sq_zeta(5) / (kappa + lambda)^2 + sq_zeta(6) / (2*kappa + lambda)^2- gamma; 
 end
 
 function df = RationalFunctionDerivative(kappa, sq_zeta, lambda)
  df = -4*sq_zeta(4)/(-2*kappa + lambda)^3 - 2*sq_zeta(5) / (kappa + lambda)^3 - 4*sq_zeta(6) / (2*kappa + lambda)^3; 
 end