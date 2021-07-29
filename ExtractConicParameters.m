% Extract conic parameters as:
% 
% a. Rotation matrix U = [u1, u2]
%
% b. Lengths along the x (u1) and y (u2) A, B (for instance, (x/A)^2 + (y/B)^2 = 1)
% 
% c. The center (position) pos of the conic.
%
% NOTE: The points transform in the canonical coordinates as x' = R*(x-pos),
%       In other words, the coordinates x in the frame of the data points are given as x = R*x' + pos

function [R, A, B, pos] = ExtractConicParameters(conic)
  
  % conic: The 6x1 or 1x6 vector of conic coeffcients
 
  ALMOST_ZERO = 1e-6;
  
  a = conic(1); b = conic(2); c = conic(3); d = conic(4); e = conic(5); f = conic(6);
  
  % The Delta quantity 
  Delta = b^2 - 4*a*c;
  absDelta = abs(Delta);
  
  % The conic PSD matrix
  C  = [ a   b/2 d/2;...
        b/2  c   e/2;...
        d/2  e/2  f];
        
  Q = C(1:2, 1:2);
  
  [U, E] = eig(Q);
  
    
  % Find r1 and r2 direction vectors as the eigenvectors of Q
  R = U;
  %if det(R) < 0, R(:, 2) = -R(:, 2), end; % make sure R is a rotation matrix 
  
  % 1. The position
  %
  % A. if Delta ~= 0 things (ellipse/hyperbola) are easy: pos = -inv(Q)*[d/2; e/2]
  v = [d/2; e/2];
  r1 = R(:, 1); 
  r2 = R(:, 2);
  lambda1 = E(1, 1);
  lambda2 = E(2, 2);
  if absDelta > ALMOST_ZERO
    pos = -inv(Q)*v;
    % To find the lengths , we solve the quadratic(s): 
    %
    % λ1*A^2 + 2*(λ1*pos'*r1+v'*r1)*A + f+v'*pos = 0 (ellipse - hyperbola A)
    % λ2*B^2 + 2*(λ2*pos'*r2+v'*r2)*B + f+v'*pos = 0 (ellipse B)
    % λ2*B^2 + 2*(λ2*pos + v)'*r2 + (sqrt(2)*r1 + r2)'*v)*B + pos'*v + 2*sqrt(2)*A*(λ1*pos + v)'*r1 + f = 0 (hyperbola B)
    % and choose the positive solutions
    %
    alpha1 = lambda1;
    beta1 = 2*(lambda1*(pos'*r1) + (v'*r1));
    gamma1 = f + v'*pos;
    D1 = beta1^2 - 4*alpha1*gamma1;
    A = abs(0.5*(-beta1 + sqrt(D1))/alpha1);
    
    if Delta < 0
      alpha2 = lambda2;
      beta2 = 2*(lambda2*(pos'*r2) + (v'*r2));
      gamma2 = f + v'*pos;
    else
      alpha2 = lambda2;  
      beta2 = 2*(lambda2*pos + v)'*r2;
      gamma2 = pos'*v + 2*sqrt(2)*A*(lambda1*pos + v)'*r1 + f;
    end
    
    D2 = beta2^2 - 4*alpha2*gamma2;
    B = abs(0.5*(-beta2 + sqrt(D2))/alpha2);
      
  else % B. Delta == 0, it's a parabola
    % NOTE: The parametric equation of the parabola is:
    %          x = A*theta^2 and y = B*theta   
    %     or,  x = A*theta   and y = B*theta^2
    %   
    if abs(lambda1) > ALMOST_ZERO
      pos = [ -v'*r1 / lambda1;...
              ((v'*r1)^2/lambda1-f)/(2*v'*r2)];
      B = abs(lambda1);
      A = sqrt(2*abs(v'*r2)); 
    else
      pos = [ ((v'*r2)^2/lambda2-f)/(2*v'*r1);...
              -v'*r2 / lambda2 ];
      B = sqrt(2*abs(v'*r1));
      A =  abs(lambda2) ; 
    end
    pos = R*pos;
  end
  
end