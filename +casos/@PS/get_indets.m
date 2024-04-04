function [x,L] = get_indets(p,I)
% Return indeterminate variables of polynomial.
%
% Available syntax:
%
%   x = get_indets(p)
%
% Returns indeterminate variables of polynomial `p`.
%
%   [...,L] = get_indets(p)
%
% If p is an n-by-m vector of polynomials with l monomial terms, this 
% returns an (n*m)-by-l logical matrix of component L_ij, indicating 
% whether the entry `p(i)` has the indeterminate variable x(j)`.
%
%   x = get_indets(p,I)
%
% Returns indeterminate variables of the polynomial expression `p(I)`.

indets = p.indets;

% get logical map of degrees
[degmat,Ldegmat] = get_degmat(p);

% map variable appearance to entries
L = logical(Ldegmat*(degmat > 0));

if nargin > 1
    % find variables that appear in subsref
    indets = indets(any(L(I,:),1));
end

x = casos.PS(indets);

end
