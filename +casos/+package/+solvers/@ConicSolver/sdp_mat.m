function m = sdp_mat(~,V,Ks,scale,dim,part)
% Index-based triangular de-vectorization for semi-definite matrices.
%
% This function takes a matrix
%
%       | v11 ... v1q |
%   V = |  :       :  |
%       | vp1 ... vpq |
%
% where each block entry vij is either a N*(N+1)/2-by-1 column or row 
% vector corresponding to stacking the lower (resp. upper) triangular
% elements of an N-by-N matrix Mij column-wise, and computes the matrix
%
%       | m11 ... m1q |
%   M = |  :       :  |
%       | mp1 ... mpq |
%
% where each block mij is a (column/row) vector that satisfies 
% mij = Mij(:).
%
% Syntax:
%
%   M = sdp_mat(V,Ks,scale,dim)
%
% Returns the block matrix M as described above with the following
% parameters:
%
% - dim:    Whether the blocks of M and V are treated as columns (dim = 1) 
%           or row vectors (dim = 2); default: dim = 1 if V is a matrix.
% - scale:  Number to de-scale the off-diagonal terms of each Mij by; 
%           default: scale = sqrt(2).
% - Ks:     Dimensions Nij of the matrices Mij; if dim = 1, then Ks is a
%           p-by-1 vector satisfying Nij = K(i); otherwise, Ks is a q-by-1
%           vector satisfying Nij = K(j) for all (i,j) in {1...p}x{1...q}.
% - part:   Optional string, either 'lower' (default) or 'upper', which
%           determines whether the lower- or upper-triangular
%           de-vectorization is performed
%

if nargin > 4
    % dimension provided, nothing to do
elseif isrow(V)
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

if nargin < 6 || isempty(part)
    part = 'lower'; % default triangular part
end
isLower = strcmpi(part,'lower'); % logical flag

% ensure matrix dimensions are a row vector
s = reshape(Ks,1,[]);

% number of elements in vectorization
Nv = s.*(s+1)/2;

assert(size(V,dim) == sum(Nv), 'Number and size of block entries in V must correspond to dimensions Ks.')

% get nonzero subindices w.r.t. block matrix V
subI = cell(1,2);
Ivec = find(sparsity(V));
[subI{:}] = ind2sub(size(V),Ivec);

% select dimension and ensure linear indices are a row vector
I = reshape(subI{dim},1,[]);

% number of elements before vectorization vij
Sv = cumsum([0 Nv(1:end-1)]);

% compute off-dimension index of corresponding block
J = sum(I > Sv', 1); % interface is 1-based

% compute linear indices of elements in each block
I0 = I - Sv(J);

% compute indices (k,l) of elements in block Vij
% ------------------------------------------------
% Convention depends on isLower:
% - If isLower = true: 
%   vectorization follows the lower-triangular part (k >= l).
%   Then I0 is given by
%       I0 = (l-1)*(2N - l + 2)/2 + (k - l + 1).
% - If isLower = false: 
%   vectorization follows the upper-triangular part (k =< l).
%   Then I0 is given by
%       I0 = (k-1)*(2N - k + 2)/2 + (l - k + 1).
%
% Note: the number of elements in the lower (resp. upper) triangle of an 
% N-by-N matrix is N*(N+1)/2.
if isLower
    l = ceil(1/2*(2*s(J) + 3 - sqrt(4*s(J).^2 + 4*s(J) + (1 - 8*I0)))) - 1;
    
    % number of rows below/including diagonal
    kprime = I0 - (l-1).*(s(J)-l/2+1);
    k = kprime + l - 1;

    % strictly lower part
    tri = (k > l);
else
    k = ceil(0.5*(2*s(J) + 3 - sqrt((2*s(J)+3).^2 - 8*I0)) - 1);

    % number of columns to the right/including diagonal
    lprime = I0 - (k-1).*(s(J)-k/2+1);
    l = lprime + k - 1;

    % strictly upper part
    tri = (k < l);
end

% de-scale strictly lower triangle
scaling = ones(size(tri));
scaling(tri) = scale;

% nonzero elements of vectorization, scaled
val = scaling.\reshape(V(Ivec),1,length(Ivec));

if nargout > 1
    % return subindices -- not implemented
    error('Not supported.')
else
    % reconstruct strictly upper triangle
    k2 = [k l(tri)]; % row
    l2 = [l k(tri)]; % col
    J2 = [J J(tri)]; % index of Mij
    Ivec2 = [Ivec Ivec(tri)]; % linear index into V

    % compute linear indices into each matrix Mij, where Im0 = N*l + k
    Im0 = s(J2).*(l2-1) + k2;
    
    % compute cumulative linear indices
    Nq = s.^2; S = cumsum([0 Nq(1:end-1)]);
    Im = S(J2) + Im0;
    
    % size of block matrix M
    sz = size(V); sz(dim) = sum(Ks.^2);
    
    % subindices and linear indices into block matrix M
    subIm{dim} = Im; subIm{3-dim} = [subI{3-dim} subI{3-dim}(tri)];
    
    % sort to match ascending indices
    [~,ii] = sort(sub2ind(sz,subIm{:}));
    
    % sparsity pattern of block matrix
    Sp = casadi.Sparsity.triplet(sz(1),sz(2),subIm{1}-1,subIm{2}-1);
    
    % return block matrix M
    m = feval(class(V),Sp,val(Ivec2(ii)));
end

end
