function [v,i,j,k,l] = sdp_vec(~,M,Ks,scale,dim,upper)
% Index-based triangular vectorization for semi-definite matrices.
%
% This function takes a matrix
%
%       | m11 ... m1q |
%   M = |  :       :  |
%       | mp1 ... mpq |
%
% where each block entry mij is either a column or row vector satisfying
% mij = Mij(:) for some N-by-N matrix Mij, and computes the matrix
%
%       | v11 ... v1q |
%   V = |  :       :  |
%       | vp1 ... vpq |
%
% where each block entry vij is a N*(N+1)/2-by-1 (column/row) vector that
% corresponds stacking the triangular elements of Mij column-wise. By
% default, the lower-triangular elements are stacked, but this can be
% optionally switched to the upper-triangular elements using the 'upper' 
% flag. 
%
% Syntax #1:
%
%   V = sdp_vec(M,Ks,scale,dim,[upper])
%
% Returns the block matrix V as described above with the following
% parameters:
%
% - dim:    Whether the blocks of M and V are treated as columns (dim = 1) 
%           or row vectors (dim = 2); default: dim = 1 if M is a matrix.
% - scale:  Number to scale the off-diagonal terms of each Mij by; 
%           default: scale = sqrt(2).
% - Ks:     Dimensions Nij of the matrices Mij; if dim = 1, then Ks is a
%           p-by-1 vector satisfying Nij = K(i); otherwise, Ks is a q-by-1
%           vector satisfying Nij = K(j) for all (i,j) in {1...p}x{1...q}.
% - upper:  Optional Boolean flag (default false), which
%           determines whether the lower- or upper-triangular elements of
%           each block Mij are stacked.
%
% Syntax #2:
%
%   [val,i,j,k,l] = sdp_vec(M,...)
%
% Returns the nonzero elements of V along with vectors of indices (i,j) and
% (k,l) corresponding to the lower-triangular elements Mij(k,l). Parameters
% apply as above.

if nargin > 4 && ~isempty(dim)
    % dimension provided, nothing to do
elseif isrow(M)
    % row vector blocks only
    dim = 2;
else
    % column-only or default for matrix
    dim = 1;
end

if nargin < 4 || isempty(scale)
    % default scaling for SCS
    scale = sqrt(2);
end

lower = nargin < 5 || ~upper; % use tril by default

% ensure matrix dimensions are a row vector
s = reshape(Ks,1,[]);
% number of elements in each matrix block
Nq = s.^2;

assert(size(M,dim) == sum(Nq), 'Number and size of block entries in M must correspond to dimensions Ks.')

% get nonzero subindices w.r.t. block matrix M
subI = cell(1,2);
[subI{:}] = ind2sub(size(M),find(sparsity(M)));

% select dimension and ensure linear indices are a row vector
I = reshape(subI{dim},1,[]);

% number of elements before matrix variables Mij
S = cumsum([0 Nq(1:end-1)]);

% compute off-dimension index of corresponding matrix variable
J = sum(I > S', 1); % interface is 1-based

% compute linear indices of elements in each matrix
I0 = I - S(J);

% compute indices (k,l) of elements in matrix Mij, where I0 = N*l + k;
% as (remainder,quotient) of linear indices divided by matrix size
l = ceil(I0 ./ s(J)); % col
k = I0 - s(J).*(l-1); % row

% determine strictly lower or upper triangle
tril_mask = (k > l); % elements below the main diagonal (row > column)
triu_mask = (k < l); % elements above the main diagonal (row < column)

% check the flag
if lower
    tri2rem = triu_mask;    % indices of strictly upper elements to remove
    tri2vec = tril_mask;    % indices of strictly lower elements to scale (off-diagonal)
else
    tri2rem = tril_mask;    % indices of strictly lower elements to remove
    tri2vec = triu_mask;    % indices of strictly upper elements to scale (off-diagonal)
end

% remove indices
J(tri2rem) = [];
l(tri2rem) = [];
k(tri2rem) = [];

% linear indices for triangular entries to keep
subItri = {subI{1}(~tri2rem) subI{2}(~tri2rem)};
Itri = sub2ind(size(M),subItri{:});

% keep tri2vec consistent with removed elements
tri2vec = tri2vec(~tri2rem);

% scale strictly the triangular elements
scaling = ones(size(tri2vec));
scaling(tri2vec) = scale;

% nonzero elements of lower triangular matrices, scaled
val = scaling.*reshape(M(Itri),1,length(Itri));

if nargout > 1
    % return subindices (i,j) of matrices Mij
    subij{dim} = J; subij{3-dim} = subItri{3-dim};
    [i,j] = subij{:};
    % return nonzero elements
    v = val;
else
    % compute linear indices into each block Vij
    % where Iv0 = (l-1)*(2N-l+2)/2 + k' and k' = k - l + 1
    if lower
        kprime = k - l + 1;                 % number of rows from diagonal
        Iv0 = (l-1).*(s(J)-l/2+1) + kprime; % number of lower-triangular elements
    else
        kprime = l - k + 1;                 % count from diagonal
        Iv0 = (k-1).*(s(J)-k/2+1) + kprime; % number of upper-triangular elements
    end

    % compute cumulative linear indices
    Nv = s.*(s+1)/2; 
    Sv = cumsum([0 Nv(1:end-1)]);
    Iv = Sv(J) + Iv0;
    
    % subindices into block matrix V
    subIv{dim} = Iv-1; 
    subIv{3-dim} = subItri{3-dim}-1; % CasADi interface has 0-index
    
    % size of block matrix V
    sz = size(M); 
    sz(dim) = sum(Nv);

    % sparsity pattern of block matrix
    Sp = casadi.Sparsity.triplet(sz(1),sz(2),subIv{:});
    
    % return block matrix V
    v = feval(class(M),Sp,val');
end

end
