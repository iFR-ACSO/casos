function [v,i,j,k,l] = sdp_vec_upper(obj, M, Ks, scale, dim)
    % Index-based upper-triangular vectorization for semi-definite matrices.

    if nargin > 3
        % dimension provided, nothing to do
    elseif isrow(M)
        % row vector blocks only
        dim = 2;
    else
        % column-only or default for matrix
        dim = 1;
    end

    if nargin < 3 || isempty(scale)
        % default scaling for SCS
        scale = sqrt(2);
    end

    % ensure matrix dimensions are a row vector
    s = reshape(Ks,1,[]);
    % number of elements in each matrix block
    Nq = s.^2;

    assert(size(M,dim) == sum(Nq), ...
        'Number and size of block entries in M must correspond to dimensions Ks.')

    % get nonzero subindices w.r.t. block matrix M
    subI = cell(1,2);
    [subI{:}]=M.sparsity().get_triplet(); % only get the indexes of the nonzero entries
    subI = {subI{1} + 1, subI{2} + 1};

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

    % determine strictly upper triangle
    triu = (k < l);

    % remove indices for strictly upper triangle
    J(triu) = [];
    l(triu) = [];
    k(triu) = [];

    % linear indices for lower triangular entries
    subItril = {subI{1}(~triu) subI{2}(~triu)};
    Itril = sub2ind(size(M),subItril{:});

    % determine strictly lower and upper triangle
    tril = (k > l);

    % scale strictly lower triangle
    scaling = ones(size(tril));
    scaling(tril) = scale;

    % nonzero elements of lower triangular matrices, scaled
    val = scaling.*reshape(M(Itril),1,length(Itril));

    if nargout > 1
        % return subindices (i,j) of matrices Mij
        subij{dim} = J; subij{3-dim} = subItril{3-dim};
        [i,j] = subij{:};
        % return nonzero elements
        v = val;
    else
        % Compute linear indices into each block Vij
        kprime = k - l + 1;  % Number of rows from diagonal
        Iv0 = (l-1).*(s(J)-l/2+1) + kprime;  % Number of lower-triangular elements

        % Compute cumulative linear indices
        Nv = s.*(s+1)/2; 
        Sv = cumsum([0 Nv(1:end-1)]);
        Iv = Sv(J) + Iv0;

        % Subindices into block matrix V
        subIv{dim} = Iv;  % 1-based indexing in MATLAB
        subIv{3-dim} = subItril{3-dim};  

        % Size of block matrix V
        sz = size(M); 
        sz(dim) = sum(Nv);
        n_rows = sz(1);
        n_cols = sz(2);

        % Convert to 0-based indices for CasADi MX triplet constructor
        row = subIv{1} - 1;
        col = subIv{2} - 1;
        
        % Create sparsity pattern from triplets
        sp = casadi.Sparsity.triplet(n_rows, n_cols, row, col);
    
        % Create sparse MX with that sparsity pattern
        V = casadi.MX(sp, val(:));

        % Return block matrix V
        v = V;
    end
end
