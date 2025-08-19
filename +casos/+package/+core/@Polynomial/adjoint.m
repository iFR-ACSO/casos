function b = adjoint(a)
% Return adjoint operator.

b = a.new_poly;

% adjoint coefficient matrix
[S,b.coeffs] = coeff_adjoint(a.get_sparsity,a.coeffs);

% set sparsity
b = set_sparsity(b,S);

end