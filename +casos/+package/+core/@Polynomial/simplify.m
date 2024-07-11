function b = simplify(a)
% Simplify polynomial expressions.

b = a.new_poly;

% simplify coefficients
coeffs = simplify(a.coeffs);

% remove zero terms
[S,b.coeffs] = coeff_update(a.get_sparsity,coeffs);

% set sparsity
b = set_sparsity(b,S);

end
