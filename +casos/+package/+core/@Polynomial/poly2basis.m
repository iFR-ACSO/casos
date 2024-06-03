function [c,S] = poly2basis(obj,S)
% Return a vector of nonzero coordinates for a given basis (sparsity).
% If no basis is given, the polynomial's sparsity pattern is used.

if nargin < 2
    % return nonzero coefficients (below)

elseif isscalar(obj) && ~isscalar(S)
    % repeat scalar inputs before projection
    obj = project(repmat(obj,size(S)),S);

else
    % project onto given sparsity pattern.
    obj = project(obj,S);
end

% copy sparsity pattern
S = sparsity(obj);

% select nonzero coefficients
c = reshape(obj.coeffs(coeff_find(S)),nnz(obj),1);

end
