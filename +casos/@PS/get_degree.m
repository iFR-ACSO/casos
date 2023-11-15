function [degrees,L] = get_degree(p,I)
% Return vector of degrees of polynomial terms.
%
% Available syntax:
%
%   degrees = get_degree(p)
%
% Returns ordered vector of degrees in polynomial `p`.
%
%   [...,L] = get_degree(p)
%
% If p is an n-by-m vector of polynomials with l monomial terms, this 
% returns an (n*m)-by-l logical matrix of component L_ij, indicating 
% whether the entry `p(i)` has term of degree `degrees(j)`.
%
%   degrees = get_degree(p,I)
%
% Returns ordered vector of degrees in the polynomial expression `p(I)`.

degsum = sum(p.degmat,2);

% vector of degrees
degrees = full(unique(degsum'));

% get logical map of degrees
[~,Ldegmat] = get_degmat(p);

% map degrees to entries
L = logical(Ldegmat*(degsum == degrees));

if nargin > 1
    % find degrees that appear in subsref
    degrees = degrees(any(L(I,:),1));
end

end
