function [S,coeffs] = coeff_blkcat(S1,S2,S3,S4,cf1,cf2,cf3,cf4)
% Block concatenate polynomial coefficient matrices.

[M,N] = size(S1);
[~,O] = size(S2);
[L,~] = size(S3);
[~,~] = size(S4);

% concatenate to a (M+L)x(NxO) matrix [S1 S2; S3 S4]
assert(size(S2,1) == M && size(S4,1) == L, 'First dimension of polynomials being concatenated are not consistent.')
assert(size(S3,2) == N && size(S4,2) == O, 'Second dimension of polynomials being concatenated are not consistent.')

% size of output
sz = [M+L, N+O];

% combine variables
dg = cell(1,4);
[indets,dg{:}] = combineVars(S1,S2,S3,S4);

% extend coefficient matrices
% let a1,...,aL and b1,...,bM be the columns of A and B
[r1,c1] = get_triplet(S1.coeffs);
[r2,c2] = get_triplet(S2.coeffs);
[r3,c3] = get_triplet(S3.coeffs);
[r4,c4] = get_triplet(S4.coeffs);

% block concatenation
% vec([A B; C D]) = [(a1 c1) ... (aN cN) ... (b1 d1) ... (bO dO)]

% horizontal offset of i-th polynomial's 1st column
off1 =     L*floor(c1/M);
off3 = M + M*floor(c3/L);
off2 =     L*floor(c2/M)  + (M+L).*N;
off4 = M + M*floor(c4/L)  + (M+L).*N;

% sparsity patterns of extended coefficient matrices
S1_coeffs = casadi.Sparsity.triplet(S1.nterm,prod(sz),r1,off1+c1);
S2_coeffs = casadi.Sparsity.triplet(S2.nterm,prod(sz),r2,off2+c2);
S3_coeffs = casadi.Sparsity.triplet(S3.nterm,prod(sz),r3,off3+c3);
S4_coeffs = casadi.Sparsity.triplet(S4.nterm,prod(sz),r4,off4+c4);

if isa(cf1,'casadi.Sparsity')
    % concatenate sparsity
    coeffs = vertcat(S1_coeffs,S2_coeffs,S3_coeffs,S4_coeffs);

else
    % extend coefficient matrices
    cf1 = sparsity_cast(cf1,S1_coeffs);
    cf2 = sparsity_cast(cf2,S2_coeffs);
    cf3 = sparsity_cast(cf3,S3_coeffs);
    cf4 = sparsity_cast(cf4,S4_coeffs);
    % concatenate coefficients
    coeffs = vertcat(cf1,cf2,cf3,cf4);
end

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, vertcat(dg{:}));

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,sz);

end
