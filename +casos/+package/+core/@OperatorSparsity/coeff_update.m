function [S,coeffs] = coeff_update(obj,coeffs,sz,dim) %#ok<INUSD>
% Update operator coefficient matrix.

if nargin > 2
    % finish after matrix operation
    error('Not supported.')
end

assert(all(size(coeffs) == [nnz(obj.sparsity_out) nnz(obj.sparsity_in)]),'Notify the developers.')

% TODO: remove zero terms
[coeffs,Si,So] = removeZero(coeffs,obj.sparsity_in,obj.sparsity_out);

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
