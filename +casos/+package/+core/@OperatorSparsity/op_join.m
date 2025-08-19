function [S,I1,I2] = op_join(obj,S2)
% Join two operator sparsity patterns.
%
% Returns indices of nonzeros in S1 and S2 relative to nonzeros in S.

assert(isequal(size(obj),size(S2)), 'Notify the developers.')

% expand coefficient matrices
[cfa,cfb,Si,So] = expand_internal(obj,S2,obj.sparsity_M,S2.sparsity_M);

% remove zero terms from addition of cfa and cfb
[coeffs,Si,So] = removeZero(cfa+cfb,Si,So);

% return sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

if nargout > 1
    % return indices
    I1 = ismembc2(find(cfa),find(coeffs));
    I2 = ismembc2(find(cfb),find(coeffs));
end

end
