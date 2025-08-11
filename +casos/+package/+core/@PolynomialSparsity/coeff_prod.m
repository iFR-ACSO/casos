function [S,coeffs] = coeff_prod(obj,coeffs,dim)
% Compute product of polynomial coefficient matrix.

nt = obj.nterm;
nv = obj.nvars;
ne = obj.numel;

% size of result
sz = sizeofMatrixOp(obj,dim);

if nt == 1
    % polynomial is a constant matrix 
    % or prod(a) = prod(c)*(x^d)^L
    if isa(coeffs,'casadi.Sparsity')
        % compute sparsity pattern
        coeffs = sparsity_prod(coeffs,size(obj),dim);

    else
        % reshape input coefficient matrix
        coeffs = prepareMatrixOp(obj,coeffs,dim);
    
        % multiply entries of coefficient
        cfb = evalMatrixOp(coeffs,@casadi_prod,dim);
    
        % reshape output coefficient matrix
        coeffs = finishMatrixOp(cfb,sz,dim);
    end

    % number of entries to multiply
    L = ne/prod(sz);
    % multiply degrees (if any)
    degmat = L*obj.degmat;

elseif is_monom(obj)
    % polynomial is a matrix of monomials
    [ia,ja] = coeff_triplet(obj);
    % match degrees to nonzero coefficients
    dga = sparse(ne,nv);
    dga(ja+1,:) = obj.degmat(ia+1,:); % NOTE: Casadi has zero index

        % combine degrees
    switch (dim)
        case 'all'
            % matrix operation over all elements
            % sum all degrees
            degmat = sum(dga,1);

        case 1
            % matrix operation along first dimension
            dgb = reshape(dga,size(obj,1),size(obj,2)*nv);
            % sum degrees columnwise
            dgc = sum(dgb,1);
            % reshape to degree matrix
            degmat = reshape(dgc,size(obj,2),nv);

        case 2
            % matrix operation along second dimension
            dgb = reshape(dga',size(obj,1)*nv,size(obj,2));
            % sum degrees rowwise
            dgc = sum(dgb,2);
            % reshape to degree matrix
            degmat = reshape(dgc,nv,size(obj,1))';

        otherwise
            error('Invalid dimension input')
    end

    % combine coefficients
    coeff1 = sum1(coeffs);

    if isa(coeffs,'casadi.Sparsity')
        % compute sparsity pattern
        cfb = sparsity_prod(coeff1,size(obj),dim);
        % diagonalize
        [~,jj] = get_triplet(cfb);
        coeffs = casadi.Sparsity.triplet(prod(sz),prod(sz),jj,jj);

    else
        % reshape input coefficient matrix
        coeffs = prepareMatrixOp(casos.Sparsity(size(obj)),coeff1,dim);
    
        % multiply coefficients
        cfb = evalMatrixOp(coeffs,@casadi_prod,dim);
    
        % expand to coefficient matrix
        coeffs = diag(cfb);
    end

    % make degree matrix unique
    [coeffs,degmat] = uniqueDeg(coeffs, degmat);

else
    % reshape input coefficient matrix
    cfa = prepareMatrixOp(obj,coeffs,dim);
    
    % number of entries to multiply
    L = ne/prod(sz);
    
    % multiply entries
    [coeffs,degmat] = prod_internal(cfa,obj.degmat,sz,dim,L,prod(sz));
end

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,obj.indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,sz);

end

function [coeffs,degmat] = prod_internal(coeffs,degmat,sz,dim,L,nb)
% Compute product of matrix elements.

nt = size(degmat,1);

if L <= 1
    % reshape output coefficient matrix
    coeffs = finishMatrixOp(coeffs,sz,dim);

    % nothing else to do
    return

else %if L == 2 || (L*nt^L) > 1e4
    % two elements to multiply or
    % number of permutations exceeds memory
    L1 = ceil(L/2);
    L2 = L - L1;
    % divide and conquer
    switch (dim)
        case {2 'all'}
            % compute product along second dimension
            cfs = blocksplit(coeffs,[0 nb*nt],[0 L1 L]);

        case 1
            % compute product along first dimension
            cfs = blocksplit(coeffs,[0 L1 L],[0 nb*nt]);

        otherwise
            error('Invalid dimension input.')
    end
    % reshape split
    cfs = [cfs{:}];
    % compute partial products
    [cfa,dga] = prod_internal(cfs{1},degmat,sz,dim,L1,nb);
    [cfb,dgb] = prod_internal(cfs{2},degmat,sz,dim,L2,nb);

    nta = size(dga,1);
    ntb = size(dgb,1);

    % multiply partial products
    % see TIMES for details
    coeffs = times_coeffs(kron(cfa,ones(ntb,1)), kron(ones(nta,1),cfb));
    degmat = times_degmat(kron(dga,ones(ntb,1)), kron(ones(nta,1),dgb));
end

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat] = removeZero(coeffs,degmat);

end
