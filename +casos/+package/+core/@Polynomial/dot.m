function c = dot(a,b)
% Dot product.

if is_operator(b)
    % composition of operators
    if ~is_operator(a), a = dualize(a); end

    assert(all(size(a.sparsity_in) == size(b.sparsity_out)), 'Dimension mismatch for operator composition.')

elseif ~is_operator(a)
    % inner product of polynomials
    assert(all(size(a) == size(b)), 'Dimension mismatch for polynomial inner product.')

else
    % evaluation of operator on polynomial
    assert(all(size(a.sparsity_in) == size(b)), 'Dimension mismatch for operator evaluation.')
end

c = a.new_poly;

% compute dot product
[S,c.coeffs] = coeff_dot(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
