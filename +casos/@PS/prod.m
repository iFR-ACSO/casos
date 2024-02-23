function B = prod(A,dim)
% Return product of array elements.

if nargin > 1
    % nothing to do
elseif isempty(A) && ~isvector(A)
    % product of empty matrix (return scalar)
    dim = 'all';
elseif isrow(A)
    % product along second dimension (return column)
    dim = 2;
else
    % product along first dimension (return row)
    dim = 1;
end

if isempty(A)
    % product of empty matrix is one
    [~,sz] = prepareMatrixOp([],A.matdim,dim);

    B = casos.PS.ones(sz);
    return

elseif isequal(dim,'all') && isscalar(A) ...
        || ~isequal(dim,'all') && size(A,dim) == 1
    % nothing to do
    B = A;
    return
end

% else:
B = casos.PS;

nt = A.nterm;
nv = A.nvars;
ne = A.numel;
sz = A.matdim;

% check number of coefficients per matrix component
idx = find(sparsity(A.coeffs));
[ia,ja] = ind2sub(size(A.coeffs),idx);

if nt == 1 
    % polynomial is a constant matrix 
    % or prod(a) = prod(c)*(x^d)^L
    
    % reshape input coefficient matrix
    [cfa,sz] = prepareMatrixOp(A.coeffs,sz,dim);

    % multiply entries of coefficient
    cfb = evalMatrixOp(cfa,@sx_prod,dim);

    % reshape output coefficient matrix
    coeffs = finishMatrixOp(cfb,sz,dim);

    % number of entries to multiply
    L = ne/prod(sz);
    % multiply degrees (if any)
    degmat = L*A.degmat;

elseif issorted(ja,'strictascend')
    % polynomial is a matrix of monomials
    % match degrees to nonzero coefficients
    dga = sparse(ne,nv);
    dga(ja,:) = A.degmat(ia,:);

    % combine degrees
    switch (dim)
        case 'all'
            % matrix operation over all elements
            % sum all degrees
            degmat = sum(dga,1);

        case 1
            % matrix operation along first dimension
            dgb = reshape(dga,sz(1),sz(2)*nv);
            % sum degrees columnwise
            dgc = sum(dgb,1);
            % reshape to degree matrix
            degmat = reshape(dgc,sz(2),nv);

        case 2
            % matrix operation along second dimension
            dgb = reshape(dga',sz(1)*nv,sz(2));
            % sum degrees rowwise
            dgc = sum(dgb,2);
            % reshape to degree matrix
            degmat = reshape(dgc,nv,sz(1))';

        otherwise
            error('Invalid dimension input')
    end

    % combine coefficients
    coeff = sum(A.coeffs,1);

    % reshape input coefficient matrix
    [cfa,sz] = prepareMatrixOp(coeff,sz,dim);

    % multiply coefficients
    cfb = evalMatrixOp(cfa,@sx_prod,dim);

    % expand to coefficient matrix
    coeffs = diag(cfb);

    % make degree matrix unique
    [coeffs,degmat] = uniqueDeg(coeffs, degmat);

else

% reshape input coefficient matrix
[cfa,sz] = prepareMatrixOp(A.coeffs,sz,dim);

% number of entries to multiply
L = ne/prod(sz);

% multiply entries
[coeffs,degmat] = prod_internal(cfa,A.degmat,sz,dim,L,prod(sz));

end

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,A.indets);

% new polynomial
B.coeffs = coeffs;
B.degmat = degmat;
B.indets = indets;
B.matdim = sz;

end

function [coeffs,degmat] = prod_internal(coeffs,degmat,sz,dim,L,nb)
% Compute product of matrix elements.

nt = size(degmat,1);

if L <= 1
    % reshape output coefficient matrix
    coeffs = finishMatrixOp(coeffs,sz,dim);

    % nothing else to do
    return

elseif L == 2 || (L*nt^L) > 1e4
    % two elements to multiply or
    % number of permutations exceeds memory
    L1 = ceil(L/2);
    L2 = L - L1;
    % divide and conquer
    switch (dim)
        case {2 'all'}
            % compute product along second dimension
            cfs = mat2cell(coeffs,nb*nt,[L1 L2]);

        case 1
            % compute product along first dimension
            cfs = mat2cell(coeffs,[L1 L2],nb*nt);

        otherwise
            error('Invalid dimension input.')
    end
    % compute partial products
    [cfa,dga] = prod_internal(cfs{1},degmat,sz,dim,L1,nb);
    [cfb,dgb] = prod_internal(cfs{2},degmat,sz,dim,L2,nb);

    nta = size(dga,1);
    ntb = size(dgb,1);

    % multiply partial products
    % see TIMES for details
    coeffs = kron(cfa,ones(ntb,1)) .* kron(ones(nta,1),cfb);
    degmat = kron(dga,ones(ntb,1)) +  kron(ones(nta,1),dgb);

else
    % (sum_a c_a[1]*x^a).* ... .* (sum_a c_a[L]*x^a)
    %  = sum_a1 ... sum_aL (c_a1[1].* ... .* c_aL[L])*x^(a1+...+aL)
    I = cell(1,L);
    [I{:}] = ndgrid(1:nt); I = cellfun(@(c) c(:),I,'UniformOutput',false);
    idx = horzcat(I{:});

    nv = size(degmat,2);
    
    % permute degree matrix
    %
    %       | dax ... day |      | da1x ... da1y        |
    %   D = |  :       :  |  ->  |   :        :     ... |
    %       | dzx ... dzy |      | daLx ... daLy        |
    %
    dga = reshape(degmat(idx',:),L,nv*nt^L);
    % combine degrees
    dgb = sum(dga,1);
    
    % permute and combine coefficients matrix
    % for details on cfa, see PREPAREMATRIXOP
    switch (dim)
        case 'all'
            % compute total product
            % rows are elements of coefficients c_a
            row = idx;              % select coefficient
            col = repmat(1:L,nt^L,1);   % select element
            ind = sub2ind(size(coeffs),row,col);
            C = reshape(coeffs(ind),nt^L,L);
    
            % combine coefficients
            cfb = sx_prod(C,2);
    
        case 2
            % compute product along second dimension
            % rows are rows of c_a, sorted by term
            row = kron((0:nb-1)',nt*ones(nt^L,L)) + repmat(idx,nb,1);
            col = repmat(1:L,nb*nt^L,1);    % select row
            ind = sub2ind(size(coeffs),row,col);
            C = reshape(coeffs(ind),nb*nt^L,L);
    
            % combine coefficients
            cfb = sx_prod(C,2);
    
        case 1
            % compute product along first dimension
            % columns are colums of c_a, sorted by coefficient
            term0 = nb*(0:nt-1)';
            row = repmat(1:L,nb*nt^L,1);    % select column
            col = kron(term0(idx),ones(nb,1)) + repmat((1:nb)',nt^L,L);
            ind = sub2ind(size(coeffs),row,col);
            C = reshape(coeffs(ind'),L,nb*nt^L);
    
            % combine coefficients
            cfb = sx_prod(C,1);

        otherwise
            error('Invalid dimension input.')
    end
    
    % reshape output degree matrix
    degmat = reshape(dgb,nt^L,nv);

    % reshape output coefficient matrix
    coeffs = finishMatrixOp(cfb,sz,dim);
end

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat] = removeZero(coeffs,degmat);

end
