% bisection routine
function [m, accuracy, iterations] = Bisection(f_handle, a, b, tolerance, maxIterations) 

% INPUT
% f : Function implementation (used to evaluate the function)
% a : left initial interval bound
% b : Right initial interval bound
% tolerance: Accuracy by which the root is approximayed
% maxIterations: Maximum number of iterations.

m = (a + b) / 2;
fa = f_handle(a);

% evaluate the function at m
fval = f_handle(m);
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
    fval = f_handle(m);
    i = i + 1;
end
iterations = i; % return number of steps
accuracy = abs(fval);
if (accuracy > tolerance)
    disp('Bisection did not converge below tolerance...')
end