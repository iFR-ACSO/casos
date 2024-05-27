function c = kron(a,b)
% Compute Kronecker product of two polynomial matrices.

% handle simple case(s) for speed up
if isempty(a) || isempty(b)
    % product with empty polynomial is empty

    c = a.new_poly(size(a).*size(b));
    return
end

% else
c = a.new_poly;

% Kronecker product of coefficients
[S,c.coeffs] = coeff_kron(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
