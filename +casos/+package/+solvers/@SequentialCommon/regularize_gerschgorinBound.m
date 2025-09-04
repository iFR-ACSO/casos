%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularizes a given Hessian matrix to ensure it is positive definite.
% The implemented method is inspired by the Levenberg method described in
% Section 2.5 in "Betts, J.- Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming, 2010". A small modification is
% done by interpolating tau parameter.
%
% Inputs:
%   H   - (n x n) Hessian matrix 
%
% Output:
%   H_reg - Regularized  Hessian matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H_reg =  regularize_gerschgorinBound(H)

epsilon = 1e-6;

% Diagonal elements of H
h_diag = diag(H);

% Sum of absolute values of elements in each row
row_sum = sum(abs(H), 2);

% Gershgorin bound
Gershgorin_bound = h_diag - (row_sum - abs(h_diag));

% get the minimum Gerschgorin bound
sigma = min(Gershgorin_bound);

% Determine tau based on sigma
if sigma > 0
    tau = 0; % Already positive definite
else
    tau = (epsilon - sigma) / (abs(sigma) + 1);
end


% Add the regularization term to the Hessian
I = eye(size(H,1)); 

% equation (2.39)
H_reg = H + tau * (abs(sigma) + 1) * I;

end
