% ========================================================================
%
%
% Test Name: unitTest_plus.m
%
% Test Description: Check if plus() returns the same coefficients
%
% Test Procedure:   Generate random coefficients for two polyomials, add
%                   them and check if coefficients are the same.
%                
% Date: 08/18/2023 
%
%
% ========================================================================

function [TestBool] = unitTest_plus()
clear
clc

% define degree and size
m   = 3;
n   = 1;
deg = 1:2;
deg2 = 2;

% generate casos monomial vector
x_cas       = casos.PS('x',m,n);
x_monom1_cas = monomials(x_cas,deg);

% generate monomial vector ordered
% [~,x_monom1_cas,~] = poly2basis(x_monom1_cas);

x_monom2_cas = monomials(x_cas,deg2);

% generate monomial vector ordered
% [~,x_monom2_cas,~] = poly2basis(x_monom2_cas);

% generate random coefficients
coeffs1 = randn(length(x_monom1_cas),1)';
coeffs2 = randn(length(x_monom2_cas),1)';

% generate casos polynomials
poly1_cas = coeffs1 * x_monom1_cas;
poly2_cas = coeffs2 * x_monom2_cas;

% adding
poly3_cas = plus(poly1_cas,poly2_cas);

% generate multipoly monomial vector as reference
x_mp         = mpvar('x',m,n);
x_monom1_mp  = monomials(x_mp,deg);

% generate sorted cell array
x_monom2_mp  = monomials(x_mp,deg2);

poly1_mp = coeffs1 * x_monom1_mp;
poly2_mp = coeffs2 * x_monom2_mp;

poly3_mp = plus(poly1_mp,poly2_mp);

% extract coefficients from new polynomial
[Vcas,~,~] = poly2basis(poly3_cas );
[Vmp,~,~]  = poly2basis(poly3_mp,monomials(x_mp,1:2));

% get double values for coefficients
poly3_cas_coeff = full(casadi.DM(Vcas));
poly3_mp_coeff  = full(Vmp);

% verify
if all(poly3_cas_coeff-poly3_mp_coeff == 0)
    TestBool = true;
else
    TestBool= false;
end
end
