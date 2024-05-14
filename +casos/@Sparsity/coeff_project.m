function [S,coeffs] = coeff_project(obj,coeffs,S)
% Project polynomial coefficient matrix onto sparsity pattern.

assert(~isa(coeffs,'casadi.Sparsity'),'Notify the developers.')

% expand coefficients to match degrees
[S,S_coeffs,coeffs] = coeff_expand(S,obj,S.coeffs,coeffs);

% project onto coefficient sparsity pattern
coeffs = project(coeffs,S_coeffs);

% remove zero terms
[coeffs,S.degmat,S.indets] = removeZero(coeffs,S.degmat,S.indets);

% store coefficients
S = set_coefficients(S,coeffs);

end
