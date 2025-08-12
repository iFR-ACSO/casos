function b = nabla(a,x)
% Evaluate polynomial nabla operator (Jacobian matrix).

assert(~is_operator(a), 'Not allowed for operators.')
assert(is_indet(x), 'Second argument must be vector of indeterminates.')

b = a.new_poly;

% compute coefficient matrix of Jacobian
[S,b.coeffs] = coeff_nabla(a.get_sparsity,a.coeffs,x);

% set sparsity
b = set_sparsity(b,S);

end
