function polyList = genRefData()

m    = 3;
n    = 1;
deg  = 1:2;
deg2 = 2;

% generate casos monomial vector
x_cas       = mpvar('x',m,n);
x_monom1_cas = monomials(x_cas,deg);

x_monom2_cas = monomials(x_cas,deg2);

% generate random coefficients
coeffs1 = randn(length(x_monom1_cas),1)';
coeffs2 = randn(length(x_monom2_cas),1)';

% generate casos polynomials
poly1_cas = coeffs1 * x_monom1_cas;
poly2_cas = coeffs2 * x_monom2_cas;


polyList = [poly1_cas;poly2_cas];
end