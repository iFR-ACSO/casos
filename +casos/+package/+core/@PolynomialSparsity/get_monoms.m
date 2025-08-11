function [z,L] = get_monoms(S,varargin)
% Return vector of monomials of sparsity pattern.
%
% Available syntax:
%
%   monoms = get_monoms(S)
%
% Returns (unsorted) vector of monomials in pattern `S`.
%
%   [...,L] = get_monoms(S)
%
% If S is an n-by-m sparsity pattern with l monomial terms, this 
% returns an (n*m)-by-l logical matrix of component L_ij, indicating 
% whether the entry `S(i)` has monomials `monoms(j)`.
%
%   monoms = get_monoms(S,I)
%
% Returns vector of monomials in the polynomial expression `p(I)`,
% if `S` is the pattern of `p`.

indets = S.indets;

% get degrees in subsref
[degmat,L] = get_degmat(S,varargin{:});

% select variables in subsref
Ivar = any(degmat > 0,1);

% monomials are already sorted canonically
z = build_monomials(degmat(:,Ivar),indets(Ivar));

end
