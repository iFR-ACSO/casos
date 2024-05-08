function S = predict_mtimes(S1,S2)
% Return sparsity pattern of matrix multiplication.

if isscalar(S1) || isscalar(S2)
    % scalar multiplication
    S = predict_times(S1,S2);
    return

elseif ~check_sz_mtimes(S1,S2)
    % dimensions are compatible if size(a,2) == size(b,1)
    throw(casos.package.core.IncompatibleSizesError.matrix(S1,S2));
end

% else:
% Repeat, permute, and reorder the coefficient matrices
%
%        | a11 ... aM1 ... a1L ... aML |
%   S1 = |  :       :       :       :  |
%        | c11 ... cM1 ... c1L ... cML |
%
% and
%
%        | x11 ... xL1 ... x1N ... xLN |
%   S2 = |  :       :       :       :  |
%        | z11 ... zL1 ... z1N ... zLN |
%
% into the element-wise product
%
%   | A |    | X |
%   | : |    | : |
%   | A |    | Z |
%   | : | .* | : |
%   | C |    | X |
%   | : |    | : |
%   | C |    | Z |
%
% where A.*X represents
%
%   | a11 ... a1L |    | x11 ... xL1 |
%   |  :       :  |    |  :       :  |
%   | aM1 ... aML |    | x11 ... xL1 |
%   |  :       :  | .* |  :       :  |
%   | a11 ... a1L |    | x1N ... xLN |
%   |  :       :  |    |  :       :  |
%   | aM1 ... aML |    | x1N ... xLN |
%

nz1 = nnz(S1);
nz2 = nnz(S2);

nt1 = S1.nterm;
nt2 = S2.nterm;

% multiply MxL with LxN matrix
[M,L] = size(S1);
[~,N] = size(S2);

% get linear indices into coefficients
% NOTE: Casadi has zero-based index
[it1,ind1] = get_triplet(S1.coeffs);
[it2,ind2] = get_triplet(S2.coeffs);

% permute rows / columns
Ind1 = repmat(ind1,1,N);
Ind2 = repmat(ind2,1,M);

% get subindices into coefficients
% NOTE: Matlab has one-based index
[i1,j1] = ind2sub([M L],Ind1+1);
[i2,j2] = ind2sub([L N],Ind2+1);

% reorder permutations
% rows of A: [a1 ... aM ... a1 ... aM]
i1 = i1 + kron(0:(N-1),M*ones(1,nz1));
% columns of B: [x1 ... x1 ... xN ... xN]
j2 = M*(j2-1) + kron(0:(M-1),ones(1,nz2)) + 1;

% permute terms
I1 = repmat(i1,1,nt2); J1 = repmat(j1,1,nt2);
I2 = repmat(i2,1,nt1); J2 = repmat(j2,1,nt1);

It1 = repmat(it1,1,N*nt2);
It2 = repmat(it2,1,M*nt1);

% reorder terms
% terms of A: [A ... A ... C ... C]
It1 = (nt2)*It1 + kron(0:(nt2-1),ones(1,N*nz1));
% terms of B: [X ... Z ... X ... Z]
It2 = It2 + kron(0:(nt1-1),nt2*ones(1,M*nz2));

% adjust terms
I1 = I1 + It1*(M*N);
J2 = J2 + It2*(M*N);

% prepare element-wise product
% NOTE: Casadi has zero-based index
cf1 = casadi.Sparsity.triplet(M*N*nt1*nt2,L,I1-1,J1-1);
cf2 = casadi.Sparsity.triplet(M*N*nt1*nt2,L,J2-1,I2-1);

% Obtain the coefficient matrix
%
%       | vec(A*X) |
%       |     :    |
%       | vec(A*Z) |
%       |     :    |
%       | vec(C*X) |
%       |     :    |
%       | vec(C*Z) |
%
% where each column vec(A*X) of S corresponds to
%
%         | <a(1,:),x(:,1)> ... <a(1,:),x(:,N)> | 
%   A*X = |        :                   :        |
%         | <a(M,:),x(:,1)> ... <a(M,:),x(:,N)> |
%
idx = find(sum2(cf1*cf2));
% transpose and reshape
[cols,rows] = ind2sub([M*N nt1*nt2],idx);
% NOTE: Casadi has zero-based index
coeffs = casadi.Sparsity.triplet(nt1*nt2,M*N,rows-1,cols-1);
% coeffs = reshape(sum2(cf1*cf2),nt1*nt2,M*N);

% return predicted sparsity pattern
S = coeff_mtimes(S1,S2,coeffs);

end
