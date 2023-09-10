function [z,L] = get_monoms(p,varargin)
% Return vector of monomials of polynomial.
%
% Available syntax:
%
%   monoms = get_monoms(p)
%
% Returns (unsorted) vector of monomials in polynomial `p`.
%
%   [...,L] = get_monoms(p)
%
% If p is an n-by-m vector of polynomials with l monomial terms, this 
% returns an (n*m)-by-l logical matrix of component L_ij, indicating 
% whether the entry `p(i)` has monomials `monoms(j)`.
%
%   monoms = get_monoms(p,I)
%
% Returns vector of monomials in the polynomial expression `p(I)`.

indets = p.indets;

% get degrees in subsref
[degmat,L] = get_degmat(p,varargin{:});

% select variables in subsref
Ivar = any(degmat > 0);

% monomials are already sorted canonically
z = build_monomials(degmat(:,Ivar),indets(Ivar));

end
