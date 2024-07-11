function [S,coeffs] = coeff_prod(S,coeffs,dim)
% Compute product of polynomial coefficient matrix.

nt = S.nterm;
nv = S.nvars;
ne = S.numel;

% size of result
sz = sizeofMatrixOp(S,dim);

if nt == 1
    % polynomial is a constant matrix 
    % or prod(a) = prod(c)*(x^d)^L
    if isa(coeffs,'casadi.Sparsity')
        % compute sparsity pattern
        coeffs = sparsity_prod(coeffs,size(S),dim);

    else
        % reshape input coefficient matrix
        coeffs = prepareMatrixOp(S,coeffs,dim);
    
        % multiply entries of coefficient
        cfb = evalMatrixOp(coeffs,@casadi_prod,dim);
    
        % reshape output coefficient matrix
        coeffs = finishMatrixOp(cfb,sz,dim);
    end

    % number of entries to multiply
    L = ne/prod(sz);
    % multiply degrees (if any)
    degmat = L*S.degmat;

elseif is_monom(S)
    % polynomial is a matrix of monomials
    [ia,ja] = coeff_triplet(S);
    % match degrees to nonzero coefficients
    dga = sparse(ne,nv);
    dga(ja+1,:) = S.degmat(ia+1,:); % NOTE: Casadi has zero index

        % combine degrees
    switch (dim)
        case 'all'
            % matrix operation over all elements
            % sum all degrees
            degmat = sum(dga,1);

        case 1
            % matrix operation along first dimension
            dgb = reshape(dga,size(S,1),size(S,2)*nv);
            % sum degrees columnwise
            dgc = sum(dgb,1);
            % reshape to degree matrix
            degmat = reshape(dgc,size(S,2),nv);

        case 2
            % matrix operation along second dimension
            dgb = reshape(dga',size(S,1)*nv,size(S,2));
            % sum degrees rowwise
            dgc = sum(dgb,2);
            % reshape to degree matrix
            degmat = reshape(dgc,nv,size(S,1))';

        otherwise
            error('Invalid dimension input')
    end

    % combine coefficients
    coeff1 = sum1(coeffs);

    if isa(coeffs,'casadi.Sparsity')
        % compute sparsity pattern
        cfb = sparsity_prod(coeff1,size(S),dim);
        % diagonalize
        [~,jj] = get_triplet(cfb);
        coeffs = casadi.Sparsity.triplet(prod(sz),prod(sz),jj,jj);

    else
        % reshape input coefficient matrix
        coeffs = prepareMatrixOp(casos.Sparsity(size(S)),coeff1,dim);
    
        % multiply coefficients
        cfb = evalMatrixOp(coeffs,@casadi_prod,dim);
    
        % expand to coefficient matrix
        coeffs = diag(cfb);
    end

    % make degree matrix unique
    [coeffs,degmat] = uniqueDeg(coeffs, degmat);

else
    % reshape input coefficient matrix
    cfa = prepareMatrixOp(S,coeffs,dim);
    
    % number of entries to multiply
    L = ne/prod(sz);
    
    % multiply entries
    [coeffs,degmat] = prod_internal(cfa,S.degmat,sz,dim,L,prod(sz));
end

% remove zero terms
[coeffs,S.degmat,S.indets] = removeZero(coeffs,degmat,S.indets);

% set dimensions of output
S.matdim = sz;
% store coefficients
S = set_coefficients(S,coeffs);

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

%% Alternative, in-situ implementation (memory-expensive)
% else
%     % (sum_a c_a[1]*x^a).* ... .* (sum_a c_a[L]*x^a)
%     %  = sum_a1 ... sum_aL (c_a1[1].* ... .* c_aL[L])*x^(a1+...+aL)
%     I = cell(1,L);
%     [I{:}] = ndgrid(1:nt); I = cellfun(@(c) c(:),I,'UniformOutput',false);
%     idx = horzcat(I{:});
% 
%     nv = size(degmat,2);
% 
%     % permute degree matrix
%     %
%     %       | dax ... day |      | da1x ... da1y        |
%     %   D = |  :       :  |  ->  |   :        :     ... |
%     %       | dzx ... dzy |      | daLx ... daLy        |
%     %
%     dga = reshape(degmat(idx',:),L,nv*nt^L);
%     % combine degrees
%     dgb = sum(dga,1);
% 
%     % permute and combine coefficients matrix
%     % for details on cfa, see PREPAREMATRIXOP
%     switch (dim)
%         case 'all'
%             % compute total product
%             % rows are elements of coefficients c_a
%             row = idx;              % select coefficient
%             col = repmat(1:L,nt^L,1);   % select element
%             ind = sub2ind(size(coeffs),row,col);
%             C = reshape(coeffs(ind),nt^L,L);
% 
%             % combine coefficients
%             cfb = sx_prod(C,2);
% 
%         case 2
%             % compute product along second dimension
%             % rows are rows of c_a, sorted by term
%             row = kron((0:nb-1)',nt*ones(nt^L,L)) + repmat(idx,nb,1);
%             col = repmat(1:L,nb*nt^L,1);    % select row
%             ind = sub2ind(size(coeffs),row,col);
%             C = reshape(coeffs(ind),nb*nt^L,L);
% 
%             % combine coefficients
%             cfb = sx_prod(C,2);
% 
%         case 1
%             % compute product along first dimension
%             % columns are colums of c_a, sorted by coefficient
%             term0 = nb*(0:nt-1)';
%             row = repmat(1:L,nb*nt^L,1);    % select column
%             col = kron(term0(idx),ones(nb,1)) + repmat((1:nb)',nt^L,L);
%             ind = sub2ind(size(coeffs),row,col);
%             C = reshape(coeffs(ind'),L,nb*nt^L);
% 
%             % combine coefficients
%             cfb = sx_prod(C,1);
% 
%         otherwise
%             error('Invalid dimension input.')
%     end
% 
%     % reshape output degree matrix
%     degmat = reshape(dgb,nt^L,nv);
% 
%     % reshape output coefficient matrix
%     coeffs = finishMatrixOp(cfb,sz,dim);
% end
% 
%%
% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat] = removeZero(coeffs,degmat);

end
