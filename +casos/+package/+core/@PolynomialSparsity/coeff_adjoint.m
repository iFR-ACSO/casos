function [S,coeffs] = coeff_adjoint(obj,coeffs)
% Return coefficient matrix of adjoint.

S = dualize(obj);

% return coefficient matrix for dualization
coeffs = sparsity_cast(coeffs,coeff_sparsity(S));

end
