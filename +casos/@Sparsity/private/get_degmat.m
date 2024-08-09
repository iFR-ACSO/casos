function [degmat,L] = get_degmat(S,I)
% Return degree matrix of sparsity pattern.
%
% Available syntax:
%
%   degmat = get_degmat(S)
%
% Returns degree matrix in pattern `S`.
%
%   [...,L] = get_degmat(S)
%
% If S is an n-by-m sparsity pattern with l monomial terms, this 
% returns an (n*m)-by-l logical matrix of component L_ij, indicating 
% whether the degree matrix of `S(i)` has a row `degmat(j,:)`.
%
%   degmat = get_degmat(S,I)
%
% Returns degree matrix in the polynomial expression `p(I)`, 
% if `S` is the pattern of `p`.

degmat = S.degmat;

% convert coefficient sparsity into logical matrix
[i,j] = coeff_triplet(S);
L = sparse(j+1,i+1,true,numel(S),S.nterm);

if nargin > 1
    % find degrees that appear in subsref
    degmat = degmat(any(L(I,:),1),:);
end

end
