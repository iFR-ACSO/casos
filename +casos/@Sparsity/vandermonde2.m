function [V,V0] = vandermonde2(S,pts,I)
% Return Vandermonde matrix for nonscalar sparsity pattern.

% L_ij is true iff S(i) has the j-th monomial
[~,L] = get_degmat(S);

% compute Vandermonde matrix for all monomials
%
%       | p1^A1 ... p1^An |
%  V' = |   :         :   |
%       | pM^A1 ... pM^An |
%
V0 = vandermonde(monomials(S),pts);

if nargin > 2
    % return Vandermonde matrix for given index
    %
    %       | p1^a1 ... p1^aN |
    %   V = |   :         :   |
    %       | pM^a1 ... pM^aN |
    %
    % where (a1,...,aN) are the monomials in S(i)
    V = V0(:,find(L(I,:))); %#ok<FNDSB>

else
    % return multidimensional Vandermonde block matrix
    %
    %       | p1^a1 ... p1^aL     0   ...   0       0   ...   0   |
    %       |   :         :       :         :       :         :   |
    %       | pL^a1 ... pM^aL     0   ...   0       0   ...   0   |
    %       |   0   ...   0     p1^b1 ... p1^bM     0   ...   0   |
    %   V = |   :         :       :         :       :         :   |
    %       |   0   ...   0     pM^b1 ... pM^bM     0   ...   0   |
    %       |   0   ...   0       0   ...   0     p1^c1 ... p1^cN |
    %       |   :         :       :         :       :         :   |
    %       |   0   ...   0       0   ...   0     pN^c1 ... pN^cN |
    %
    % where (a1,...,aL), (b1,...,bM), (c1,...,cN) are the monomials in each
    % element of S, respectively.
    ne = S.numel;
    % elements in V0
    nV = numel(V0);

    % number of terms in each element
    nts = sum(L,2);
    % offset of each element's Vandermonde block
    off = cumsum([0; nts(1:end-1)]);

    % select block matrix (all terms)
    I = kron((1:ne)',ones(nV,1));
    % iterate nonzero entries of V (all terms)
    J = repmat((1:nV)',ne,1);
    % convert to row and column subindices
    [i0,j0] = ind2sub(size(V0),J);

    % dimension of each block corresponds to number of terms
    tf = (i0 <= nts(I)) & (j0 <= nts(I));
    % remove additional terms and points
    I(~tf) = [];

    % row indices of nonzero entries of V
    ii = off(I) + i0(tf);
    % column indices of nonzero entries of V
    jj = off(I) + j0(tf);

    % repeat monomial assignment matrix
    LL = kron(L',true(length(pts),1));
    % find monomials of each element
    [jL,iL] = find(LL);
    % convert to row and column subindices
    [i1,j1] = ind2sub(size(V0),jL);

    % number of rows corresponds to number of evaluation points
    tf = (i1 <= nts(iL));
    % remove additional points
    i1(~tf) = [];
    j1(~tf) = [];
    % select elements in total Vandermonde matrix
    idx = sub2ind(size(V0),i1,j1);

    % build block matrix
    V = casadi.DM.triplet(ii-1,jj-1,V0(idx),sum(nts),sum(nts));
end

end
