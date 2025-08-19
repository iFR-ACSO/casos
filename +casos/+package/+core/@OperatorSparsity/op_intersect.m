function [S,I1,I2] = op_intersect(obj,S2)
% Intersect two operator sparsity patterns.
%
% Returns indices of nonzeros in S relative to nonzeros in S1 and S2.

assert(isequal(size(obj),size(S2)), 'Notify the developers.')

% expand coefficient matrices
[cfa,cfb,Si,So] = expand_internal(obj,S2,obj.sparsity_M,S2.sparsity_M);

% intersection of cfa and cfb
coeffs = cfa*cfb;

if nargout > 1
    % return indices
    I1 = ismembc2(find(coeffs),find(cfa));
    I2 = ismembc2(find(coeffs),find(cfb));
end

% remove zero terms from intersection
[coeffs,Si,So] = removeZero(coeffs,Si,So);

% return sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
