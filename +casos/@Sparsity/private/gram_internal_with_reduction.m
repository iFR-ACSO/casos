function [Z,K,Mp,Md] = gram_internal_with_reduction(Lz,degmat,indets, Lz_red)
% Compute Gram basis and mappings.

lp = size(Lz,1);

% dimension K(i) of Gram basis for p(i)
K = full(sum(Lz_red,2));

% build square matrix of monomials
nt = size(degmat,1);
% vectorized square matrix S = z*z^T
D = kron(degmat,ones(nt,1)) + kron(ones(nt,1),degmat);
% build logical map for square matrix
% -> L(i,j) is true iff square matrix for p(i) includes S(j)
% that is, S(j) = z(j1)*z(j2) 
% and Gram form of p(i) includes z(j1) and z(j2)
L = kron(Lz,ones(1,nt)) & kron(ones(1,nt),Lz);
Lr = kron(Lz_red,ones(1,nt)) & kron(ones(1,nt),Lz_red);

% remove unused monomials from square matrix
I = any(L,1);
L(:,~I) = [];
Lr(:,~I) = [];
D(~I,:) = [];

% build coefficients for output
nT = size(D,1);

% enumerate entries of S, first for p(1), then p(2), ...
% index into entries -> set {(i,j) : p(j) includes S(i)}
[i,j] = find(L');

% coefficient matrix for Z
coeffs = casadi.Sparsity.triplet(nT,lp,i-1,j-1);

% set output
Z = casos.Sparsity;
[Z.coeffs,Z.degmat,~] = uniqueDeg(coeffs,D);
Z.indets = indets;
Z.matdim = [lp 1];

% for the reduced version
[i,j] = find(Lr');
coeffs = casadi.Sparsity.triplet(nT,lp,i-1,j-1);
Zr = casos.Sparsity;
[Zr.coeffs,Zr.degmat,Iuni] = uniqueDeg(coeffs,D);
Zr.indets = indets;
Zr.matdim = [lp 1];

% compute primal mapping
[~,ic] = Iuni{:};

% enumerate nonzero elements in L' relative to nonzeros of Luni'
idx = ismembc2(sub2ind(size(Z.coeffs),ic(i),j), find(Z.coeffs));

% create mapping from nonzero elements of Lr to nonzero elements Lruni
Mp = sparse(idx,1:nnz(Lr),1,nnz(Z.coeffs),nnz(Lr));

% return adjoint inverse
row_sums = sum(Mp, 2);
row_sums(row_sums == 0) = 1; 
Md = spdiags(1./row_sums, 0, size(Mp,1), size(Mp,1)) * Mp;
 
end



