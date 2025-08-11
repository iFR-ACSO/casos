function r = coeff_properint(obj,coeffs)
% Compute integral over [0 1] as double.

% integral of monomials
mprod = 1./prod(obj.degmat+1,2);

% multiply with coefficients + reshape
r = full(reshape(mprod'*coeffs, size(obj)));

end
