function s = pnorm2(a)
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

s = properintc(a'*a);

end