% ========================================================================
%
%
% Test Name: unitTest_monomials.m
%
% Test Description: Generate monomial vector
%
% Test Procedure:   A monomial vector is generated using casos and compare 
%                   to monomial vector generated by mulitpoly toolbox.
%                
% Date: 08/18/2023 
%
%
% ========================================================================

function [TestBool] = unitTest_monomials()

% define degree and size
m   = 3;
n   = 1;
deg = 1:2;

% generate casos monomial vector
x_cas       = casos.PS('x',m,n);
x_monom_cas = monomials(x_cas,deg);

% generate sorted cell array
x_monom_cas_str = str(x_monom_cas);

% generate multipoly monomial vector as reference
x_mp        = mpvar('x',m,n);
x_monom_mp  = monomials(x_mp,deg);

% generate sorted cell array
x_monom_mp_str = char(x_monom_mp);

% verify
if all(strcmp(x_monom_cas_str,x_monom_mp_str))
    TestBool = true;
else
    TestBool = false;
end

end