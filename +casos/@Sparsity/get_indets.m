function [x,L] = get_indets(S,I)
% Return indeterminate variables of sparsity pattern.
%
% Available syntax:
%
%   x = get_indets(S)
%
% Returns indeterminate variables of pattern `S`.
%
%   [...,L] = get_indets(S)
%
% If S is an n-by-m sparsity pattern with l monomial terms, this 
% returns an (n*m)-by-l logical matrix of component L_ij, indicating 
% whether the entry `S(i)` has the indeterminate variable x(j)`.
%
%   x = get_indets(S,I)
%
% Returns indeterminate variables of the polynomial expression `p(I)`,
% if `S` is the pattern of `p`.

indets = S.indets;

% get logical map of degrees
[degmat,Ldegmat] = get_degmat(S);

% map variable appearance to entries
L = logical(Ldegmat*(degmat > 0));

if nargin > 1
    % find variables that appear in subsref
    indets = indets(any(L(I,:),1));
end

x = casos.PS(indets);

end
