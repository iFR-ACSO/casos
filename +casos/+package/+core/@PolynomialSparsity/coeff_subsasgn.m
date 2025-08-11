function [S,coeffs] = coeff_subsasgn(obj,S2,coeff1,coeff2,ii)
% Subassignment of polynomial coefficient matrices.

% combine variables
[indets,dga,dgb] = combineVar(obj.indets,S2.indets,obj.degmat,S2.degmat);

% extend rhs coefficient matrix
[i2,j2] = coeff_triplet(S2);
S2_coeffs = casadi.Sparsity.triplet(S2.nterm,obj.numel,i2,ii(j2+1)-1); % NOTE: Casadi has zero-based index

% remove selected lhs coefficients
[i1,j1] = coeff_triplet(obj); I = ismember(j1+1,ii);
S1_coeffs = casadi.Sparsity.triplet(obj.nterm,obj.numel,i1(~I),j1(~I));

if isa(coeff1,'casadi.Sparsity')
    % assignment
    coeffs = vertcat(S1_coeffs,S2_coeffs);

else
    % extend coefficient matrices
    cf1 = project(coeff1,S1_coeffs);
    cf2 = sparsity_cast(coeff2,S2_coeffs);
    % assigment
    coeffs = vertcat(cf1,cf2);
end

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, [dga;dgb]);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,size(obj));

end
