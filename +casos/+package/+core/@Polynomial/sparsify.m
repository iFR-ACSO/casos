function b = sparsify(obj)
% Sparsify polynomial expressions.

b = obj.new_poly;

% sparsify coefficients
coeffs = sparsify(obj.coeffs);

% remove zero terms
[S,b.coeffs] = coeff_update(obj.get_sparsity,coeffs);

% set sparsity
b = set_sparsity(b,S);

end
