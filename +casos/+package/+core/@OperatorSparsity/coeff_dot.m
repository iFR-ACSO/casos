function [S,coeffs] = coeff_dot(obj,S2,coeff1,coeff2)
% Compute coefficient matrix for dot operation.

assert(is_operator(S2), 'Notify the developers.')
assert(numel(obj) == numel(S2), 'Inputs must be of compatible size.')

% expand coefficient matrices
[cf1,cf2] = expand_internal(obj,S2,coeff1,coeff2);

% return dot product of coefficients
coeffs = dot(cf1,cf2);

% dot product of polynomials is always dense scalar
S = casos.Sparsity.scalar;

end
