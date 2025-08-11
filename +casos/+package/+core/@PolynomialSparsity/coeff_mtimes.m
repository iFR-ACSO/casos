function [S,coeffs] = coeff_mtimes(obj,S2,coeffs)
% Matrix multiplication of two polynomial coefficient matrices.

nta = obj.nterm;
ntb = S2.nterm;

% combine variables
[indets,dga,dgb] = combineVar(obj.indets,S2.indets,obj.degmat,S2.degmat);

% (sum_a c_a*x^a)*(sum_b c_b*x^b) = (sum_a sum_b (c_a*c_b)*(x^a*x^b)
degmat = times_degmat(kron(dga,ones(ntb,1)), kron(ones(nta,1),dgb));

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,[size(obj,1) size(S2,2)]);

end
