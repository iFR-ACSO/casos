function [S,coeffs] = coeff_project_operator(obj,coeffs,Si,So,keep_zeros)
% Project operator coefficient matrix onto input-output patterns.

assert(~isa(coeffs,'casadi.Sparsity'),'Notify the developers.')

% join input/output nonzero coordinates
[Si2,I1i] = op_join(obj.sparsity_in,Si);
[So2,I1o] = op_join(obj.sparsity_out,So);

% expand matrix to joint nonzeros
coeff0 = expand_matrix(obj.sparsity_M,coeffs,[nnz(So2) nnz(Si2)],I1i,I1o);

% find common nonzero coordinates
[~,I2i,~] = op_intersect(Si2,Si);
[~,~,I2o] = op_intersect(So2,So);

% return nonzeros
coeffs = coeff0(I2o,I2i);

if nargin < 4 || ~keep_zeros
    % remove zero terms
    [coeffs,Si,So] = removeZero(coeffs,Si,So);
end

% return operator sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
