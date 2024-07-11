function b = transpose(a)
% Transpose of polynomial matrix.

b = a.new_poly;

% transpose coefficient matrix
[S,b.coeffs] = coeff_transpose(a.get_sparsity,a.coeffs);

% set sparsity
b = set_sparsity(b,S);

end
