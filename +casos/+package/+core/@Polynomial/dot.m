function c = dot(a,b)
% Dot product.

c = a.new_poly;

if is_operator(a) && ~is_operator(b)
    % evaluate dual operator (linear form)
    assert(is_dual(a), 'Only allowed for duals. Use "evaluate" instead.')

    assert(all(size(a.sparsity_in) == size(b)), 'Dimension mismatch for operator evaluation.')

    % evaluation of operator on polynomial
    [S,c.coeffs] = coeff_evaluate(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

else
    % inner product of polynomials
    assert(is_operator(a) == is_operator(b), 'Must not mix polynomials and operators.')
    assert(all(size(a) == size(b)), 'Dimension mismatch for inner product.')

    % compute dot product
    [S,c.coeffs] = coeff_dot(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);
end

% set sparsity
c = set_sparsity(c,S);

end
