function c = mtimes(a,b)
% Multiply (intersect) two polynomial sparsity patterns.

a = casos.Sparsity(a);
b = casos.Sparsity(b);

% sparsity patterns must be of same size
if ~check_sz_equal(a,b)
    throw(casos.package.core.IncompatibleSizesError.other(a,b));
end

% expand coefficient matrices
[c,cfa,cfb] = coeff_expand(a,b,a.coeffs,b.coeffs);

% remove zero terms
[c.coeffs,c.degmat,c.indets] = removeZero(cfa*cfb,c.degmat,c.indets);

end
