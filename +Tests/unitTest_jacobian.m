% ========================================================================
%
%
% Test Name: unitTest_poly2jacobian.m
%
% Test Description: Check if jacobian returns the same results
%
% Test Procedure:   Generate monomial vector and a corresponding
%                   polynomial calculate jacobian. Compare coefficients to
%                   check validity.
%                
% Date: 08/18/2023 
%
%
% ========================================================================

function [TestBool] = unitTest_jacobian()
clear
clc

% define degree and size
m   = 3;
n   = 1;
deg = 1:2;

% generate casos monomial vector
x_cas       = casos.PS('x',m,n);
x_monom_cas = monomials(x_cas,deg);

% bring monomials in compareable order
[~,x_monom_cas,~] = poly2basis(x_monom_cas);

% generate casos polynomial
coeffs1 = randn(length(x_monom_cas),1)';

% generate polynomial
poly_cas = coeffs1*x_monom_cas;

jac_cas = jacobian(poly_cas,x_cas);

% generate multipoly monomial vector as reference
x_mp        = mpvar('x',m,n);
x_monom_mp  = monomials(x_mp,deg);

% generate multipoly polynomial
poly_mp = coeffs1*x_monom_mp;

jac_mp= jacobian(poly_mp,x_mp);

% verification: check each polynomial in jacobian
for k = 1:length(jac_cas)

% get coefficients
[Vcas,~,~] = poly2basis(jac_cas(k));
[Vmp,~,~]  = poly2basis(jac_mp(k),monomials(x_mp,jac_mp.mindeg:jac_mp.maxdeg));

  % check if equal value apart numerical differences
  if ~all(full(casadi.DM(Vcas))-full(Vmp) <= 1e-12) % avoid numerical issues; 10^(-12) is zero enough
    TestBool = false;
    return
  end

end

TestBool = true;

end
