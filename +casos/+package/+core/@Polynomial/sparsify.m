function b = sparsify(a)
% Sparsify polynomial expressions.

b = a.new_poly;

% sparsify coefficients
coeffs = sparsify(a.coeffs);

% remove zero terms
[S,b.coeffs] = coeff_update(a.get_sparsity,coeffs);

% set sparsity
b = set_sparsity(b,S);

end
