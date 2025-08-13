function S = new_from_coefficients(coeffs,Si,So)
% Create sparsity pattern with (generic) coefficients.

if ~isa(coeffs,'casadi.Sparsity')
    coeffs = sparsity(coeffs);
end

% new sparsity pattern
S = casos.Sparsity.create(casos.package.core.OperatorSparsity(coeffs,Si,So));

end
