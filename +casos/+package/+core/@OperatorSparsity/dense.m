function S = dense(varargin)
% Create dense operator pattern.

assert(nargin > 0, 'Not enough input arguments.')
assert(nargin < 3, 'Too many input arguments.')

if isnumeric(varargin{1})
    % input-output dimensions given
    Si = casos.Sparsity.dense(varargin{end},1);
    So = casos.Sparsity.dense(varargin{1},1);

else
    % input-output patterns given
    Si = casos.Sparsity(varargin{1});
    So = casos.Sparsity(varargin{end});
end

% coefficients are dense
M = casadi.Sparsity.dense(nnz(So),nnz(Si));

% create operator pattern
S = casos.package.core.OperatorSparsity(M,Si,So);

end
