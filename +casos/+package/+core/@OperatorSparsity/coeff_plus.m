function [S,coeffs] = coeff_plus(obj,S2,coeff1,coeff2)
% Add two operator coefficient matrices.

% expand internally
[cfa,cfb,Si,So] = expand_internal(obj,S2,coeff1,coeff2);

% remove zero terms
[coeffs,Si,So] = removeZero(cfa+cfb,Si,So);

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
