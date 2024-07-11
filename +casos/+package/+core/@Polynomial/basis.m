function S = basis(p,I)
% Return basis of multidimensional polynomial.
%
% Note: This function is made obsolete by the use of polynomial sparsity
% and only kept for use with indexing argument.

if nargin < 2
    % fall back to sparsity
    warning('Deprecated: Use sparsity() instead.')
    S = sparsity(p);

else
    % return sparsity of subindex
    S = sub(p.get_sparsity,find(I),casos.Sparsity(nnz(I),1));
end

end
