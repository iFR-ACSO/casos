function [S,coeffs] = coeff_mtimes(S1,S2,coeffs)
% Matrix multiplication of two polynomial coefficient matrices.

S = casos.Sparsity;

nta = S1.nterm;
ntb = S2.nterm;

% combine variables
[indets,dga,dgb] = combineVar(S1.indets,S2.indets,S1.degmat,S2.degmat);

% (sum_a c_a*x^a)*(sum_b c_b*x^b) = (sum_a sum_b (c_a*c_b)*(x^a*x^b)
degmat = times_degmat(kron(dga,ones(ntb,1)), kron(ones(nta,1),dgb));

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity pattern
S.degmat = degmat;
S.indets = indets;
S.matdim = [size(S1,1) size(S2,2)];
% store coefficients
S = set_coefficients(S,coeffs);

end
