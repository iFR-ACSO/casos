% ========================================================================
%
%
% Test Name: unitTest_degrees.m
%
% Test Description: Check if degrees() returns the same min. and max deg.
%
% Test Procedure:   Generate polyomial with casos and multipoly and compare
%                   results
%
%                
% Date: 08/18/2023 
%
%
% ========================================================================

function [TestBool] = unitTest_degrees()
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
poly_cas = casos.PS.sym('c',x_monom_cas,n);

[min_cas,max_cas] = degrees(poly_cas);


% degrees only for casos

% % generate multipoly monomial vector as reference
% x_mp        = mpvar('x',m,n);
% x_monom_mp  = monomials(x_mp,deg);
% 
% % generate multipoly polynomial
% poly_mp = polydecvar('c1',x_monom_mp);

% verify
if (min_cas == deg(1)) && (max_cas == deg(2))
    TestBool = true;
else
    TestBool= false;
end
end
