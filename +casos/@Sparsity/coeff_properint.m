function r = coeff_properint(S,coeffs)
% Compute integral over [0 1] as double.

% integral of monomials
mprod = 1./prod(S.degmat+1,2);

% multiply with coefficients + reshape
r = full(reshape(mprod'*coeffs, size(S)));

end
