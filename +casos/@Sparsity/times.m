function c = times(a,b)
% Return sparsity pattern of element-wise multiplication.

a = casos.Sparsity(a);
b = casos.Sparsity(b);

% sparsity patterns must be of same size
if ~check_sz_equal(a,b)
    throw(casos.package.core.IncompatibleSizesError.other(a,b));
end

% join coefficient matrices
c = coeff_times(a,b,a.coeffs,b.coeffs);

end
