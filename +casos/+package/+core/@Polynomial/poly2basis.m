function [c,S] = poly2basis(obj,varargin)
% Return a vector of nonzero coordinates for a given basis (sparsity).
% If no basis is given, the polynomial's sparsity pattern is used.
%
% Deprecated: Unless used with indexing argument.

assert(~is_operator(obj), 'Not allowed for operators.')

if nargin < 2 || ~islogical(varargin{1})
    % use coordinates
    warning('Deprecated: Use coordinates() instead.')
    [c,S] = coordinates(obj,varargin{:});
    return
end

% indexing argument
I = varargin{1};

if isempty(obj) || all(~I)
    % empty polynomial, selection, or projection
    assert(~isempty(obj) || all(~I),'Cannot project empty polynomial.')

    c = sparse(0,1);
    S = casos.Sparsity(0,0);

else
    % nonzeros and basis of subreference
    [i,j] = coeff_triplet(obj.get_sparsity);
    tf = I(j+1);
    idx = sub2ind(size(obj.coeffs),i(tf)+1,j(tf)+1);
    c = reshape(obj.coeffs(idx),nnz(tf),1);
    S = basis(obj,I);
end

end
