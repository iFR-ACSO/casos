function [Z,K,Mp,Md] = gram_internal(Lz,degmat,indets)
% Compute Gram basis and mappings.

lp = size(Lz,1);

% dimension K(i) of Gram basis for p(i)
K = full(sum(Lz,2));

% build square matrix of monomials
nt = size(degmat,1);
% vectorized square matrix S = z*z^T
D = kron(degmat,ones(nt,1)) + kron(ones(nt,1),degmat);
% build logical map for square matrix
% -> L(i,j) is true iff square matrix for p(i) includes S(j)
% that is, S(j) = z(j1)*z(j2) 
% and Gram form of p(i) includes z(j1) and z(j2)
L = kron(Lz,ones(1,nt)) & kron(ones(1,nt),Lz);

% remove unused monomials from square matrix
I = any(L,1);
L(:,~I) = [];
D(~I,:) = [];

% build coefficients for output
%
%     [ S(i1), ..., S(iN) ]
% Z = [ S(j1), ..., S(jM) ]
%     [ S(k1), ..., S(kL) ]
%
% where i_, j_, k_ are those indices satisfying that
% L(1,i_), L(2,j_), L(3,k_), respectively, are true
nT = size(D,1);

% enumerate entries of S, first for p(1), then p(2), ...
% index into entries -> set {(i,j) : p(j) includes S(i)}
[i,j] = find(L');

% coefficient matrix for Z
coeffs = casadi.Sparsity.triplet(nT,lp,i-1,j-1);

% set output
Z = casos.Sparsity;
[Z.coeffs,Z.degmat,Iuni] = uniqueDeg(coeffs,D);
Z.indets = indets;
Z.matdim = [lp 1];

% compute primal mapping
[id,ic] = Iuni{:};
% enumerate nonzero elements in L' relative to nonzeros of Luni'
% NOTE: We use an internal MATLAB function here for speed-up 
% since the result of find is already sorted.
idx = ismembc2(sub2ind(size(Z.coeffs),ic(i),j), find(Z.coeffs));
% create mapping from nonzero elements of L to nonzero elements Luni
Mp = sparse(idx,1:nnz(L),1,nnz(Z.coeffs),nnz(L));

% enumerate nonzero elements in L'
% [il,jl] = find(Luni');
% Innz = sparse(il,jl,1:nnz(Luni),length(id),lp);
% create mapping from nonzero elements of L to nonzero elements Luni
% ii = sub2ind(size(Luni'),ic(i),j);
% Mp = sparse(Innz(ii),1:nnz(L),1,nnz(Luni),nnz(L));

% return adjoint inverse
Md = sum(Mp,2).\Mp;

end
