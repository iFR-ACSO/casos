function g = nabla(f,x)
% Evaluate polynomial nabla operator (Jacobian matrix).

assert(is_indet(x), 'Second argument must be vector of indeterminates.')

g = f.new_poly;

% compute coefficient matrix of Jacobian
[S,g.coeffs] = coeff_nabla(f.get_sparsity,f.coeffs,x);

% set sparsity
g = set_sparsity(g,S);

end
