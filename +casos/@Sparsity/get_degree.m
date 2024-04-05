function [degrees,L] = get_degree(S,I)
% Return vector of degrees of sparsity pattern.
%
% Available syntax:
%
%   degrees = get_degree(S)
%
% Returns ordered vector of degrees in pattern `S`.
%
%   [...,L] = get_degree(S)
%
% If S is an n-by-m sparsity pattern with l monomial terms, this 
% returns an (n*m)-by-l logical matrix of component L_ij, indicating 
% whether the entry `S(i)` has term of degree `degrees(j)`.
%
%   degrees = get_degree(S,I)
%
% Returns ordered vector of degrees in the polynomial expression `p(I)`,
% if `S` is the pattern of `p`.

degsum = sum(S.degmat,2);

% vector of degrees
degrees = full(unique(degsum'));

% get logical map of degrees
[~,Ldegmat] = get_degmat(S);

% map degrees to entries
L = logical(Ldegmat*(degsum == degrees));

if nargin > 1
    % find degrees that appear in subsref
    degrees = degrees(any(L(I,:),1));
end

end
