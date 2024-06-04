function [c,S] = poly2basis(obj,S)
% Return a vector of nonzero coordinates for a given basis (sparsity).
% If no basis is given, the polynomial's sparsity pattern is used.

if nargin < 2
    % return nonzero coefficients (below)
    S = sparsity(obj);
    coeffs = obj.coeffs;

elseif isempty(obj) || (isscalar(obj) && isempty(S))
    % empty polynomial
    assert(isempty(S),'Cannot project empty polynomial.')

    c = sparse(0,1);
    return

elseif isscalar(obj) && ~isscalar(S)
    % repeat scalar inputs before projection
    [sp_rep,coeffs] = coeff_repmat(obj.get_sparsity,obj.coeffs,size(S));
    coeffs = coeff_project(sp_rep,coeffs,S,true);

else
    % project onto given sparsity pattern.
    [~,coeffs] = coeff_project(obj.sparsity,obj.coeffs,S,true);
end

% select nonzero coefficients
c = reshape(coeffs(coeff_find(S)),nnz(S),1);

end
