function s = pnorm2(a)
<<<<<<< HEAD
% function S = pnorm2(A)
%
% DESCRIPTION
%   Computes squared polynomial norm as multiple integral over [0 1].
%
% INPUTS
%   A: polynomial
%   X: vector of polynomial variables [default: X = A.varname]
%
% OUTPUTS
%   B: polynomial
%

% 06/10/2022 TC  Initial Coding

s = properint(a'*a);
=======
% Compute polynomial 2-norm as integral of a^2 over [0 1].
%
% Note: This 2-norm does NOT correspond to the norm induced by dot.

% TODO: use internal operations
b = (a'*a);

% integral of monomials
mprod = 1./prod(b.degmat+1,2);

% multiply with coefficients + reshape
s = full(reshape(mprod'*b.coeffs, b.matdim));
>>>>>>> origin/main

end