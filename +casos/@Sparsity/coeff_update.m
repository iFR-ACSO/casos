function [S,coeffs] = coeff_update(S,coeffs,sz,dim)
% Update polynomial coefficient matrix.

if nargin > 3
    % finish after matrix operation
    coeffs = finishMatrixOp(coeffs,sz,dim);
end

assert(all(size(coeffs) == [S.nterm prod(sz)]),'Notify the developers.')

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,S.degmat,S.indets);

% new sparsity pattern
S.degmat = degmat;
S.indets = indets;
S.matdim = sz;
% store coefficients
S = set_coefficients(S,coeffs);

end
