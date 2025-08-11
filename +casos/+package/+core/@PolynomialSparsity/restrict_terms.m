function [S,I] = restrict_terms(S,deg)
% Restrict monomials terms to specified degrees.

% choose terms to keep
I = ismember(S.degsum,deg);

% restrict coefficients
coeffs = sub(S.coeffs,find(I)-1,(1:S.numel)-1);
% restrict degree matrix
degmat = S.degmat(I,:);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,S.indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,size(S));

end
