%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The implemented method regularizes the Hessian by computung the frobenius
% norm and to ensure positive definitness. A "new" Hessian is than
% constructed. This method is based on [1] p. 49-50.
%
% [1] - Nocedal, J and Wright, S.J. - Numerical Optimization, 
%       Second Edition,  Springer Series in Operations Research
%       and Financial Engineering
%
% Inputs:
%   H   - (n x n) Hessian matrix 
%
% Output:
%   H_reg - Regularized Hessian matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H_reg = regularize_minFrobNorm(H)

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

% [1] equation (3.43)
tau = zeros(length(v),1);
idx = find(v < epsilon) ;

tau(idx) = epsilon - v(idx);

% get correction matrix
D_corr = diag(tau); 

H_reg = V * (D + D_corr) * V';   

% Ensure symmetry
H_reg = (H_reg + H_reg') / 2;

end
