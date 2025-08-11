function [S,coeffs] = coeff_plus(obj,S2,cfa,cfb)
% Add (join) two polynomial coefficient matrices.

% combine variables
[indets,dga,dgb] = combineVar(obj.indets,S2.indets,obj.degmat,S2.degmat);

% make degree matrix unique
[coeffs,degmat] = uniqueDeg([cfa;cfb], [dga;dgb]);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,size(obj));

end
