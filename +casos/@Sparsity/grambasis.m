function [Z,K,z,Mp,Md] = grambasis(S,I)
% Return Gram basis of polynomial vector.

if nargin < 2
    lp = numel(S);
    I = true(lp,1);
    idx = 1:lp;
else
    lp = nnz(I);
    idx = find(I);
end

% get logical maps for degrees and indeterminate variables
% -> Ldeg(i,j) is true iff p(i) has terms of degree(j)
% -> Lvar(i,j) is true iff p(i) has terms in indets(j)
[degree,Ldeg] = get_degree(S,I);
[indets,Lvar] = get_indets(S,I);
% remove non-indexed logicals
Ldeg(~I,:) = []; Lvar(~I,:) = [];
% detect degrees and variables
Id = any(Ldeg,1); Iv = any(Lvar,1);
% remove unused degrees or variables
Ldeg(:,~Id) = []; Lvar(:,~Iv) = [];

% min and max degree of Gram basis vector
mndg = floor(min(degree)/2);
mxdg =  ceil(max(degree)/2);

% ensure min and max degree are even
[degree,ii] = unique([degree,(2*mndg):(2*mxdg)]);
ld = length(degree); tmp = [Ldeg false(lp,ld)];
Ldeg = tmp(:,ii);

% split into even and odd monomials
Ldeg_e = Ldeg(:,mod(degree,2)==0);
Ldeg_o = Ldeg(:,mod(degree,2) >0);
% add even degrees if necessary
Ldeg_e(:,1:end-1) = Ldeg_e(:,1:end-1) | Ldeg_o;
Ldeg_e(:,2:end) = Ldeg_e(:,2:end) | Ldeg_o;

% get half-degree vector of monomials
% TODO: compute degree matrix directly (avoid unnecessary checks)
z = to_vector(casos.Sparsity.scalar(indets,mndg:mxdg));
% compute logical map for half-degree monomials
% -> Lz_deg(i,j) is true iff z(i) is of degree(j)/2
% -> Lz_var(i,j) is true iff z(i) includes indets(j)
[~,Lz_deg] = get_degree(z);
[~,Lz_var] = get_indets(z);

% build logical map for half-degree monomials
% -> Lz(i,j) is true iff Gram form of p(i) includes z(j)
% that is, z(j) is of a degree in the Gram basis for p(i)
% AND z(j) only has indeterminate variables that p(i) has, too
Lz = (Ldeg_e * Lz_deg' & ~(~Lvar * Lz_var'));

% discard monomials based on simple checks
[~,Ldegmat] = get_degmat(S);
lz = numel(z);
% perform checks vector-wise
MX = arrayfun(@(i) repmat(max(S.degmat(Ldegmat(i,:),Iv),[],1),lz,1), idx, 'UniformOutput', false);
MN = arrayfun(@(i) repmat(min(S.degmat(Ldegmat(i,:),Iv),[],1),lz,1), idx, 'UniformOutput', false);
% Imx = any(  ceil(max(p.degmat,[],1)/2) < z.degmat , 2 );
% Imn = any( floor(min(p.degmat,[],1)/2) > z.degmat , 2 );
zdm = repmat(z.degmat,lp,1);
Irem = [ceil(vertcat(MX{:})/2) < zdm, floor(vertcat(MN{:})/2) > zdm];
% Lz(:,Imx | Imn) = false;
Lz(reshape(any(Irem,2),lz,lp)') = false;

% remove unused monomials from base vector
I = any(Lz,1);
Lz(:,~I) = [];
degmat = z.degmat(I,:);
z = to_vector(build_monomials(degmat,z.indets)); % do we need this?

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
[Z.coeffs,Z.degmat,~,Mp] = uniqueDeg(coeffs,D);
Z.indets = z.indets;
Z.matdim = [lp 1];

% return adjoint inverse
Md = sum(Mp,2).\Mp;

end
