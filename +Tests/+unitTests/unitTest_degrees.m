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

% define degree and size
m   = 3;
n   = 1;
deg = 1:2;

% generate casos monomial vector
x_cas       = casos.PS('x',m,n);
x_monom_cas = monomials(x_cas,deg);

% generate casos polynomial
poly_cas = casos.PS.sym('c',x_monom_cas,n);

[min_cas,max_cas] = degrees(poly_cas);



% verify
if (min_cas == deg(1)) && (max_cas == deg(2))
    TestBool = true;
else
    TestBool= false;
end
end
