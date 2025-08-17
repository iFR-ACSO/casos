function [S,coeffs] = coeff_project(obj,coeffs,S2,keep_zeros)
% Project operator coefficient matrix onto sparsity pattern.

assert(~isa(coeffs,'casadi.Sparsity'),'Notify the developers.')

% expand coefficients to match input-output patterns
[S_coeffs,coeffs,Si,So] = expand_internal(S2,obj,S2.coeff_sparsity,coeffs);

if isa(coeffs,'casadi.DM') && ~is_regular(coeffs)
    % handle irregular coefficients (inf, nan)
    warning('Projection of irregular coefficients not yet implemented.')
end

% project onto coefficient sparsity pattern
coeffs = project(coeffs,S_coeffs);

if nargin < 4 || ~keep_zeros
    % remove zero terms
    [coeffs,Si,So] = removeZero(coeffs,Si,So);
end

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
