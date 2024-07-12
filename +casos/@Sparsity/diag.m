function S = diag(n,k,varargin)
% Create diagonal sparsity pattern.

if nargin < 2
    % diag(n)
    S = casos.Sparsity(casadi.Sparsity.diag(n));
elseif isnumeric(k)
    % diag(n,k,...)
    S = casos.Sparsity(casadi.Sparsity.diag(n,k),varargin{:});
else
    % diag(n,...)
    S = casos.Sparsity(casadi.Sparsity.diag(n),k,varargin{:});
end

end
