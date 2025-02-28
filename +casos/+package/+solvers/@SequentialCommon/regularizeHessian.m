%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularizes a given Hessian matrix to ensure it is positive definite.
% The implemented method regularizes the Hessian by computung the frobenius
% norm and to ensure positive definitness. A "new" Hessian is than
% constructed.
% 
% Inputs:
%   H   - (n x n) Hessian matrix (must be symmetric)
%   tau - Regularization parameter (small positive scalar, e.g., 1e-6)
% 
% Output:
%   H_reg - Regularized positive definite Hessian matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H_reg = regularizeHessian(H)
% regularizeHessian: Regularizes a given Hessian matrix to ensure it is positive definite.
% 
% Inputs:
%   H   - (n x n) Hessian matrix (must be symmetric)
%   tau - Regularization parameter (small positive scalar, e.g., 1e-6)
% 
% Output:
%   H_reg - Regularized positive definite Hessian matrix


% Check if the input matrix is symmetric
if ~issymmetric(H)
   H = (H + H') / 2;
end

% Perform Eigenvalue Decomposition (EVD)
[Q, D] = eig(H); % Q: Eigenvectors, D: Diagonal matrix of eigenvalues

% Extract eigenvalues from diagonal matrix D
lambda = diag(D);

epsilon = 1e-6;

% use Frobenius norm of current approximation
tau = epsilon * norm(H, 'fro');

% Regularize eigenvalues: replace negative or small values with tau
lambda_reg = max(lambda, tau);

% Reconstruct the regularized Hessian
D_reg = diag(lambda_reg); % Create the diagonal matrix of regularized eigenvalues
H_reg = Q * D_reg * Q';   % Reconstruct the Hessian using eigenvectors and regularized eigenvalues

% Ensure symmetry (numerical precision may cause small asymmetries)
H_reg = (H_reg + H_reg') / 2;

end
