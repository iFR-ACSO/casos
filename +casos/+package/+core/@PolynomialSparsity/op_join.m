function [S,I1,I2] = op_join(S1,S2)
% Join two polynomial sparsity patterns.
%
% Returns indices of nonzeros in S1 and S2 relative to nonzeros in S.

assert(isequal(size(S1),size(S2)), 'Notify the developers.')

% expand coefficient matrices
[cfa,cfb,degmat,indets] = expand_internal(S1,S2,S1.coeffs,S2.coeffs);

% remove zero terms from addition of cfa and cfb
[coeffs,degmat,indets] = removeZero(cfa+cfb,degmat,indets);

% return sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,size(S1));

if nargout > 1
    % return indices
    I1 = ismembc2(find(cfa),find(coeffs));
    I2 = ismembc2(find(cfb),find(coeffs));
end

end
