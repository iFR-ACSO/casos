function c = dot(a,b)
% Dot product.

if is_operator(b)
    % dot product of operators
    a = dualize(a);

elseif ~is_operator(a)
    % dot product of polynomials
    assert(numel(obj) == numel(S2), 'Inputs must be of compatible size.')
end

c = a.new_poly;

% compute dot product
[S,c.coeffs] = coeff_dot(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
