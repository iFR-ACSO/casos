function s = pnorm2(a)
% Compute polynomial 2-norm as integral of a^2 over [0 1].
%
% Note: This 2-norm does NOT correspond to the norm induced by dot.

% TODO: use internal operations
b = (a'*a);

% integral of monomials
mprod = 1./prod(b.degmat+1,2);

% multiply with coefficients + reshape
s = full(reshape(mprod'*b.coeffs, b.matdim));

end