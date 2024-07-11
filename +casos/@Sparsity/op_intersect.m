function [S,I1,I2] = op_intersect(S1,S2)
% Intersect two polynomial sparsity patterns.
%
% Returns indices of nonzeros in S relative to nonzeros in S1 and S2.

assert(isequal(size(S1),size(S2)), 'Notify the developers.')

S = casos.Sparsity;

% expand coefficient matrices
[cfa,cfb,degmat,indets] = expand_internal(S1,S2,S1.coeffs,S2.coeffs);

% intersection of cfa and cfb
coeffs = cfa*cfb;

if nargout > 1
    % return indices
    I1 = ismembc2(find(coeffs),find(cfa));
    I2 = ismembc2(find(coeffs),find(cfb));
end

% remove zero terms from intersection
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% return sparsity pattern
S.degmat = degmat;
S.indets = indets;
S.matdim = size(S1);
% store coefficients
S = set_coefficients(S,coeffs);

end
