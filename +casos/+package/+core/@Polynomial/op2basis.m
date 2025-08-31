function [M,S] = op2basis(obj,varargin)
% Return a matrix of nonzero coordinates for given input and output basis
% (sparsity). If no basis is given, the operator's sparsity patterns are
% used.
%
% Deprecated: Unless used with input-output sparsity patterns.

assert(is_operator(obj), 'Only allowed for operators.')

if nargin < 3
    % use coordinates
    warning('Deprecated. Use coordinates() instead.')
    [c,S] = coordinates(obj,varargin{:});
    % return matrix
    M = sparsity_cast(c,coeff_sparsity(S));
    return
end

% input-output sparsity patterns
Si = varargin{1};
So = varargin{2};

assert(isequal(size(Si),size(obj.sparsity_in)),'Input dimensions must not change.')
assert(isequal(size(So),size(obj.sparsity_out)),'Output dimensions must not change.')

% project onto input-output sparsity
[S,M] = coeff_project_operator(obj.get_sparsity,obj.coeffs,Si,So,true);

end
