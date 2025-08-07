%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularizes a given Hessian matrix to ensure it is positive definite.
% The implemented method is inspired by the Levenberg method described in
% Section 2.5 in "Betts, J.- Practical Methods for Optimal Control and 
% Estimation Using Nonlinear Programming, 2010". A small modification is 
% done by interpolating tau parameter.
% 
% Inputs:
%   H   - (n x n) Hessian matrix (must be symmetric)
%   tau - Regularization parameter (small positive scalar, e.g., 1e-6)
% 
% Output:
%   H_reg - Regularized positive definite Hessian matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H_regularized = hessian_Levenberg(H)

    epsilon = sqrt(eps);

    % Diagonal elements of H
    h_diag = diag(H);
    
    % Sum of absolute values of elements in each row
    row_sum = sum(abs(H), 2);
    
    % Gershgorin bound (vectorized)
    Gershgorin_bound = h_diag - (row_sum - abs(h_diag));

    % Take the minimum Gershgorin bound 
    sigma = min(Gershgorin_bound);
    
     % Determine tau based on sigma
    if sigma > 0
        tau = 0; % Already positive definite
    else
        tau = (epsilon - sigma) / (abs(sigma) + 1);
    end

    
    % Add the regularization term to the Hessian
    I = eye(size(H,1)); % Identity matrix
    H_regularized = H + tau * (abs(sigma) + 1) * I;
end
