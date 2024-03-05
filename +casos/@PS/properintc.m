function d = properintc(a)
% Integration in all variables from zero to one, respectively.

a = casos.PS(a);

% get polynomial info
acoef = a.coeffs;
adeg  = a.degmat;
sza   = a.matdim;

% integral of monomials
mprod = 1./prod(adeg+1,2);

% multiply with coefficients + reshape
d = full(reshape(mprod'*acoef, sza));

end
