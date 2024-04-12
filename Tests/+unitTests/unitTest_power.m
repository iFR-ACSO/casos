% ========================================================================
%
%
% Test Name: unitTest_power.m
%
% Test Description: Check if power() returns the same coefficients
%
% Test Procedure:   Generate random coefficients for two polyomials,
%                   multiply them and check if coefficients are the same.
%                
% Date: 08/18/2023 
%
%
% ========================================================================

function [TestBool] = unitTest_power()

% define degree and size
m   = 3;
n   = 1;
deg = 1:2;
deg2 = 2;

% generate casos monomial vector
x_cas       = casos.PS('x',m,n);
x_monom1_cas = monomials(x_cas,deg);


% generate random coefficients
coeffs1 = randn(length(x_monom1_cas),1)';

% generate casos polynomials
poly1_cas = coeffs1 * x_monom1_cas;

% adding
polyPow_cas = power(poly1_cas,3);

% generate multipoly monomial vector as reference
x_mp         = mpvar('x',m,n);
x_monom1_mp  = monomials(x_mp,deg);

poly1_mp = coeffs1 * x_monom1_mp;

polyPow_mp = power(poly1_mp,3);

% extract coefficients from new polynomial
[Vcas,~,~]   = poly2basis(polyPow_cas );
[Vmp,~,~]    = poly2basis(polyPow_mp,monomials(x_mp,polyPow_mp.mindeg:polyPow_mp.maxdeg));

% get double values for coefficients
poly_cas_coeff = full(casadi.DM(Vcas));
poly_mp_coeff  = full(Vmp);

% verify
if all(poly_cas_coeff-poly_mp_coeff <= 1e-12) % avoid numerical issues; 10^(-12) is zero enough
    TestBool = true;
else
    TestBool= false;
end
end
