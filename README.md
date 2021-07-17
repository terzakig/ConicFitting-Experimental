# ConicFitting-Experimental
Optimal conic fitting

An approach for solving the problem of fitting conics to sparse data without having to solve the egeneralized eigenvalue poroblem. Instead, the idea here is to use a projection operator for 6D vectors onto the ellipse/parabola/hyperbola manifold and thereafter use them as starting point for local search to optimal solution. This solution is thereafter improved to reach a geometrically unbiased solution (i.e. use nearest point-to-conic errors instead of the usual algebraic errors) 
