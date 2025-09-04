%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The implemented method regularizes the Hessian by computung the frobenius
% norm and to ensure positive definitness. A "new" Hessian is than
% constructed.
%
% Inputs:
%   H   - (n x n) Hessian matrix 
%
% Output:
%   H_reg - Regularized Hessian matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H_reg = regularize_scaledFrobNorm(H)

% Check if the input matrix is symmetric
if ~issymmetric(H)
    H = (H + H') / 2;
end

% Perform Eigenvalue Decomposition (EVD) 
[V, D] = eig(H); 

% get eigenvalues from (diagonal) matrix
v = diag(D);

% tolerance
epsilon = 1e-6;

% use Frobenius norm of current approximation
tau = max(epsilon, epsilon * norm(H, 'fro'));

% clipping
lambda_reg = max(tau,v);

% Reconstruct the "clipped" Hessian
Lambda_reg = diag(lambda_reg); 
H_reg = V * Lambda_reg * V';   

% Ensure symmetry
H_reg = (H_reg + H_reg') / 2;

end
