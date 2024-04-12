% ========================================================================
%
%
% Test Name: unitTest_poly2basis.m
%
% Test Description: Check if poly2basis returns the same results
%
% Test Procedure:   Generate monomial vector and a corresponding
%                   polynomial and check if poly2basis returns the same 
%                   monomials and coefficients.
%                
% Date: 08/18/2023 
%
%
% ========================================================================

function [TestBool] = unitTest_poly2basis()

% define degree and size
m   = 3;
n   = 1;
deg = 1:2;

% generate casos monomial vector
x_cas       = casos.PS('x',m,n);
x_monom_cas = monomials(x_cas,deg);


% generate casos polynomial
poly_cas = casos.PS.sym('c',x_monom_cas,n);

% generate multipoly monomial vector as reference
x_mp        = mpvar('x',m,n);
x_monom_mp  = monomials(x_mp,deg);

% generate multipoly polynomial
poly_mp = polydecvar('c1',x_monom_mp);

% pol2basis of casos and multipoly
[~,Rcas,~] = poly2basis(poly_cas,x_monom_cas);
[~,Rmp,~]  = poly2basis(poly_mp,x_monom_mp);

% check if same monomials are returned
Rcas_str = str(Rcas);
Rcas_mp  = char(Rmp);

% test
if all(strcmp( Rcas_str,Rcas_mp  ))
    TestBool = true;
else
    TestBool= false;
end

end
