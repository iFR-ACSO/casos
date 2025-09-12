function v = coeff_getnz(obj,coeffs,I)
% Return nonzero elements of a polynomial coefficients matrix.

if isa(coeffs,'casadi.Sparsity')
    % return dense sparsity pattern
    v = casadi.Sparsity.dense(numel(I));
    return
end

% else:
nz = sparsity_cast(coeffs,casadi.Sparsity.dense(nnz(obj),1));

% select nonzeros
v = nz(I);

end
