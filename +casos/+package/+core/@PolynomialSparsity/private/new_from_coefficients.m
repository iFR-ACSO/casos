function S = new_from_coefficients(coeffs,degmat,indets,matdim)
% Create sparsity pattern with (generic) coefficients.

if ~isa(coeffs,'casadi.Sparsity')
    coeffs = sparsity(coeffs);
end

% new sparsity pattern
S = casos.Sparsity.create(casos.package.core.PolynomialSparsity(coeffs,degmat,indets,matdim));

end
