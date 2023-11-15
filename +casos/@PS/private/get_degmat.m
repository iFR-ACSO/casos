function [degmat,L] = get_degmat(p,I)
% Return degree matrix of polynomial.
%
% Available syntax:
%
%   degmat = get_degmat(p)
%
% Returns degree matrix in polynomial `p`.
%
%   [...,L] = get_degmat(p)
%
% If p is an n-by-m vector of polynomials with l monomial terms, this 
% returns an (n*m)-by-l logical matrix of component L_ij, indicating 
% whether the degree matrix of `p(i)` has a row `degmat(j,:)`.
%
%   degmat = get_degmat(p,I)
%
% Returns degree matrix in the polynomial expression `p(I)`.

degmat = p.degmat;

% sparsity of coefficients
S = sparsity(p.coeffs);
% convert CasADi sparsity into logical matrix
[i,j] = ind2sub(size(S),find(S));
L = sparse(j,i,true,size2(S),size1(S));

if nargin > 1
    % find degrees that appear in subsref
    degmat = degmat(any(L(I,:),1),:);
end

end
