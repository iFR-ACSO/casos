function [S,coeffs] = coeff_repmat(obj,coeffs,rep)
% Repeat copies of operator.

assert(isrow(rep) && length(rep) == 2, 'Replication factors must be a pair (row) of integers or two integer scalars.')

% repeat input-output patterns
Si = repeat_sparsity(obj.sparsity_in,rep(2));
So = repeat_sparsity(obj.sparsity_out,rep(1));

% repeat coefficients
coeffs = repmat(coeffs,rep);

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end

function S = repeat_sparsity(S,n)
% Repeat polynomial sparsity pattern.

if iscolumn(S)
    % repeat column by first dimension
    S = repmat(S,[n 1]);

else
    % repeat row (or matrix) by second dimension
    S = repmat(S,[1 n]);
end

end
