function Z = basis(p,I)
% Return base matrix of multidimensional polynomial.

if nargin < 2
    I = true(size(p));
    lp = numel(p);
else
    lp = nnz(I);
end

% degree matrix and logical map for monomials
% -> L(i,j) is true iff p(i) has monomial z(j)
[D,L] = get_degmat(p,I);

L(~I,:) = [];
% remove unused monomials
Jd = any(L,1);
L(:,~Jd) = [];
% remove unused variables
Jv = any(D~=0,1);
D(:,~Jv) = [];

% build coefficients for output
%
%       [ z(i1) ... z(iN) |   0   ...   0   |   0   ...   0   | ... ]
% Z^T = [   0   ...   0   | z(j1) ... z(jM) |   0   ...   0   | ... ]
%       [   0   ...   0   |   0   ...   0   | z(k1) ... z(kL) | ... ]
%
% where i_, j_, k_ are those indices satisfying that
% L(1,i_), L(2,j_), L(3,k_), respectivel, are true
lZ = sum(L,"all");
nT = size(D,1);

% enumerate entries of z, first for p(1), then p(2), ...
% index into entries -> set {(i,j) : p(j) includes z(i)}
[i,j] = find(L');

% coefficient matrix for Z
coeffs = casadi.SX.triplet(i-1, (j(:)-1)'.*lZ+(1:lZ)-1, ones(lZ,1), nT, lp*lZ);

% set output
Z = casos.PS;
[Z.coeffs,Z.degmat] = uniqueDeg(coeffs,D);
Z.indets = p.indets(Jv);
Z.matdim = [lZ lp];
