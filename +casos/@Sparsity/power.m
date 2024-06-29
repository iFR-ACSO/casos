function b = power(a,n)
% Add (join) two polynomial sparsity patterns.

a = casos.Sparsity(a);

% sparsity pattern and exponent must be of same size
if ~check_sz_equal(a,n)
    throw(casos.package.core.IncompatibleSizesError.other(a,n));
end

% power of coefficient matrix
b = coeff_power(a,a.coeffs,n);

end
