% ========================================================================
%
%
% Test Name: unitTest_kron.m
%
% Test Description: Check if kron() returns the same coefficients
%
% Test Procedure:   Generate random coefficients for two polyomials,
%                   calculate the Kronecker product them and check 
%                   if coefficients are the same.
%                
% Date: 08/18/2023 
%
%
% ========================================================================

function [TestBool] = unitTest_kron()

% define degree and size
m   = 3;
n   = 1;
deg = 1:2;
deg2 = 2;

% generate casos monomial vector
x_cas        = casos.PS('x',m,n);
x_monom1_cas = monomials(x_cas,deg);


x_monom2_cas = monomials(x_cas,deg2);


% generate random coefficients
coeffs1 = randn(length(x_monom1_cas),1)';
coeffs2 = randn(length(x_monom2_cas),1)';

% generate casos polynomials
poly1_cas = coeffs1 * x_monom1_cas;
poly2_cas = coeffs2 * x_monom2_cas;

% adding
poly3_cas = kron(poly1_cas,poly2_cas);

% generate multipoly monomial vector as reference
x_mp         = mpvar('x',m,n);
x_monom1_mp  = monomials(x_mp,deg);

% generate sorted cell array
x_monom2_mp  = monomials(x_mp,deg2);

poly1_mp = coeffs1 * x_monom1_mp;
poly2_mp = coeffs2 * x_monom2_mp;

poly3_mp = kron(poly1_mp,poly2_mp);

% extract coefficients from new polynomial
[Vcas,~,~]   = poly2basis(poly3_cas );
[Vmp,~,~]  = poly2basis(poly3_mp,monomials(x_mp,poly3_mp.mindeg:poly3_mp.maxdeg));

% get double values for coefficients
poly3_cas_coeff = full(casadi.DM(Vcas));
poly3_mp_coeff  = full(Vmp);

% verify
if all(poly3_cas_coeff-poly3_mp_coeff <= 1e-12) % avoid numerical issues; 10^(-12) is zero enough
    TestBool = true;
else
    TestBool= false;
end
end