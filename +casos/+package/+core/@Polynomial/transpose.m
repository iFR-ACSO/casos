function b = transpose(a)
% Transpose of polynomial matrix.

assert(~is_operator(a), 'Not allowed for operators. Use "adjoint" instead.')

b = a.new_poly;

% transpose coefficient matrix
[S,b.coeffs] = coeff_transpose(a.get_sparsity,a.coeffs);

% set sparsity
b = set_sparsity(b,S);

end
