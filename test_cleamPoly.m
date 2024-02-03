clc
clear

deg = 0:2;
m   = 3;
n   = 1;

% generate casos monomial vector
x_cas       = casos.PS('x',m,n);
x_monom     = monomials(x_cas,deg);

% generate random coefficient between an lower and upper bound
L = 1e-10;
U = 1e-6;

coeffs = U + (L-U).*randn(length(x_monom),1)';

% polynomial
poly = coeffs*x_monom

Asopt = cleanpoly(to_multipoly(poly),1e-6)
A     = cleanpoly(poly,1e-6)



B     = cleanpoly(poly,[], 2:3)
Bsopt = cleanpoly(to_multipoly(poly),[], 2:3)


C = cleanpoly(poly,[],{'x_1' 1:2})
Csopt = cleanpoly(to_multipoly(poly),[],{'x_1' 1:2})
