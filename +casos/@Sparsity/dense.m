function S = dense(n,k,varargin)
% Create dense sparsity pattern.

if nargin < 2
    % dense(n)
    S = casos.Sparsity(casadi.Sparsity.dense(n));
elseif isnumeric(k)
    % dense(n,k,...)
    S = casos.Sparsity(casadi.Sparsity.dense(n,k),varargin{:});
else
    % dense(n,...)
    S = casos.Sparsity(casadi.Sparsity.dense(n),k,varargin{:});
end

end
