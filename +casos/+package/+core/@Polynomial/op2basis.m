function [M,S] = op2basis(obj,varargin)
% Return a matrix of nonzero coordinates for given input and output basis
% (sparsity). If no basis is given, the operator's sparsity patterns are
% used.
%
% Deprecated.

warning('Deprecated. Use coordinates() instead.')

assert(is_operator(obj), 'Only allowed for operators.')

% get coordinates
[c,S] = coordinates(obj,varargin{:});

% return matrix
M = sparsity_cast(c,coeff_sparsity(S));

end
