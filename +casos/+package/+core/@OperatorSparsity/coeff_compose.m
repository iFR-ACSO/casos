function [S,coeffs] = coeff_compose(obj,S2,coeff1,coeff2)
% Compute coefficient matrix for operator composition.

assert(is_operator(S2), 'Notify the developers.')
assert(numel(obj.sparsity_in) == numel(S2.sparsity_out), 'Dimensions mismatch for composition.')

% find common nonzero coordinates
[~,I1,I2] = op_intersect(obj.sparsity_in,S2.sparsity_out);

% coordinate matrices
M1 = get_submatrix(coeff1,1:nnz(obj.sparsity_out),I1);
M2 = get_submatrix(coeff2,I2,1:nnz(S2.sparsity_in));

if isa(M1,'casadi.Sparsity')
    % matrix product on sparsity patterns
    %
    % compute the matrix product
    %
    %   | a11 ... a1L |   | b11 ... b1N |
    %   |  :       :  | * |  :       :  |
    %   | aM1 ... aML |   | bL1 ... bLN |
    %
    % via the element-wise product
    %
    %   | a11 ... a1L |    | b11 ... bL1 |
    %   |  :       :  |    |  :       :  |
    %   | aM1 ... aML |    | b11 ... bL1 |
    %   |  :       :  | .* |  :       :  |
    %   | a11 ... a1L |    | b1N ... bLN |
    %   |  :       :  |    |  :       :  |
    %   | aM1 ... aML |    | b1N ... bLN |
    %   
    A = repmat(M1,size(M2,2),1);
    B = T(reshape(repmat(M2,size(M1,1),1),size(A,2),size(A,1)));

    coeffs = reshape(sum2(A*B),size(M1,1),size(M2,2));

else
    % composition of operator matrices
    coeffs = M1*M2;
end

% remove zero terms
[coeffs,Si,So] = removeZero(coeffs,S2.sparsity_in,obj.sparsity_out);

% new operator pattern
S = new_from_coefficients(coeffs,Si,So);

end
