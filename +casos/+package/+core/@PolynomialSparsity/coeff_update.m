function [S,coeffs] = coeff_update(obj,coeffs,sz,dim)
% Update polynomial coefficient matrix.

if nargin > 3
    % finish after matrix operation
    coeffs = finishMatrixOp(coeffs,sz,dim);
elseif nargin < 3
    % matrix dimensions do not change
    sz = size(obj);
end

assert(all(size(coeffs) == [obj.nterm prod(sz)]),'Notify the developers.')

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,obj.degmat,obj.indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,sz);

end
