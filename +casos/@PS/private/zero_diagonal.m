function Lz = zero_diagonal(pdegmat, degmat, Lz)
% zero diagonal algorithm (similar in practive to simplify from sosopt)

% build square matrix of monomials
nt = size(degmat,1);

% vectorized square matrix S = z*z^T
D = kron(degmat,ones(nt,1)) + kron(ones(nt,1),degmat);

% build logical map for square matrix
L = kron(Lz,ones(1,nt)) & kron(ones(1,nt),Lz);

% zz^t matrix
zzdegmat = D(L,:);

% obtain only diagonal elements
dzdegmat = zzdegmat(1:sum(Lz)+1:sum(L), :);
idx = 1:sum(Lz);

% degree of the diagonal elements must be unique and no other way to obtain
% those degrees should be possible
[~, B1] = arrayfun(@(i) ismember(dzdegmat(i,:), zzdegmat, 'rows'), idx, 'UniformOutput', false);
[~, B2] = arrayfun(@(i) ismember(dzdegmat(i,:), zzdegmat, 'rows', 'legacy'), idx, 'UniformOutput', false);
B1 = horzcat(B1{:});
B2 = horzcat(B2{:});
uni = B1 == B2;

% verify if the degree associated with a diagonal element exists in the
% degree basis of p
ediag = arrayfun(@(i) ismember(dzdegmat(i,:), pdegmat, 'rows'), idx, 'UniformOutput', false);
ediag = ~horzcat(ediag{:});

% does not exist in p and the diagonal degree is unique 
ediag = ediag & uni;

% where Lz is 1
Lz_ones = find(Lz==1);

% indiate which monomial is not necessary
Lz(Lz_ones(ediag==1)) = 0;

end

