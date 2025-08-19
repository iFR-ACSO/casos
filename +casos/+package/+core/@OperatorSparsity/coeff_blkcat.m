function [S,coeffs] = coeff_blkcat(obj,S2,S3,S4,cf1,cf2,cf3,cf4)
% Coefficient matrix of block concatenation.

assert(is_matrix(obj) && is_matrix(S2) && is_matrix(S3) && is_matrix(S4), 'Not implemented.')

% the block concatenation
%
%   | A  B |
%   | C  D |
%
% corresponds to the operator
%
%   (x,y) -> (Ax+By, Cx+Dy)
%
[M,N] = size(obj);
[~,O] = size(S2);
[L,~] = size(S3);
[~,~] = size(S4);

% concatenate to a (M+L)x(NxO) matrix [S1 S2; S3 S4]
assert(size(S2,1) == M && size(S4,1) == L, 'First dimension of polynomials being concatenated are not consistent.')
assert(size(S3,2) == N && size(S4,2) == O, 'Second dimension of polynomials being concatenated are not consistent.')

% join inputs to S1 and S3
[Si1,I1i,I3i] = op_join(obj.sparsity_in,S3.sparsity_in);
% join inputs to S2 and S4
[Si2,I2i,I4i] = op_join(S2.sparsity_in,S4.sparsity_in);
% join outputs of S1 and S2
[So1,I1o,I2o] = op_join(obj.sparsity_out,S2.sparsity_out);
% join outputs of S3 and S4
[So2,I3o,I4o] = op_join(S3.sparsity_out,S4.sparsity_out);

% expand coefficient matrices
cfA = expand_matrix(obj.sparsity_M,cf1,[nnz(So1) nnz(Si1)],I1i,I1o);
cfB = expand_matrix(S2.sparsity_M,cf2,[nnz(So1) nnz(Si2)],I2i,I2o);
cfC = expand_matrix(S3.sparsity_M,cf3,[nnz(So2) nnz(Si1)],I3i,I3o);
cfD = expand_matrix(S4.sparsity_M,cf4,[nnz(So2) nnz(Si2)],I4i,I4o);

% concatenate input/output sparsity patterns
Si = vertcat(Si1,Si2);
So = vertcat(So1,So2);

% block concatenate coefficient matrices
coeffs = blockcat(cfA,cfB,cfC,cfD);

% remove zero terms
[coeffs,Si,So] = removeZero(coeffs,Si,So);

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
