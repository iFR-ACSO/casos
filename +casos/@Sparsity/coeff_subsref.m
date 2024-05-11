function [S,coeffs] = coeff_subsref(obj,coeffs,ii,sz)
% Subreference polynomial coefficient matrix.

nt = obj.nterm;

S = casos.Sparsity;

if isa(coeffs,'casadi.Sparsity')
    % subreference sparsity pattern
    coeffs = sub(coeffs,1:nt,ii,true);

else
    % subreference coefficient matrix
    coeffs = coeffs(:,ii);
end

% remove coefficients, degrees, and/or indeterminates 
% that do not appear in the referenced polynomial
[coeffs,S.degmat,S.indets] = removeZero(coeffs,obj.degmat,obj.indets);
% resize
S.matdim = sz;
% store coefficients
S = set_coefficients(S,coeffs);

end
