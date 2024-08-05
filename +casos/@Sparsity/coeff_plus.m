function [S,coeffs] = coeff_plus(S1,S2,cfa,cfb)
% Add (join) two polynomial coefficient matrices.

S = casos.Sparsity;

% combine variables
[indets,dga,dgb] = combineVar(S1.indets,S2.indets,S1.degmat,S2.degmat);

% make degree matrix unique
[coeffs,degmat] = uniqueDeg([cfa;cfb], [dga;dgb]);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity pattern
S.degmat = degmat;
S.indets = indets;
S.matdim = size(S1);
% store coefficients
S = set_coefficients(S,coeffs);

end
