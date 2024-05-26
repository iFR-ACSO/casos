function S = dense(n,k,varargin)
% Create dense sparsity pattern.

if isnumeric(k)
    % diag(n,k,...)
    S = casos.Sparsity(casadi.Sparsity.dense(n,k),varargin{:});
else
    % diag(n,...)
    S = casos.Sparsity(casadi.Sparsity.dense(n),k,varargin{:});
end

end
