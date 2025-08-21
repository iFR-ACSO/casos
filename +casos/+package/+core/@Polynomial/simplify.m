function b = simplify(obj)
% Simplify polynomial expressions.

b = obj.new_poly;

% simplify coefficients
coeffs = simplify(obj.coeffs);

% remove zero terms
[S,b.coeffs] = coeff_update(obj.get_sparsity,coeffs);

% set sparsity
b = set_sparsity(b,S);

end
