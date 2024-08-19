function c = kron(a,b)
% Compute Kronecker product of two polynomial sparsity patterns.

a = casos.Sparsity(a);
b = casos.Sparsity(b);

% handle simple case(s) for speed up
if isempty(a) || isempty(b)
    % product with empty polynomial is empty

    c = casos.Sparsity(size(a).*size(b));
    return
end

% else
% Kronecker product of coefficients
c = coeff_kron(a,b,a.coeffs,b.coeffs);

end
