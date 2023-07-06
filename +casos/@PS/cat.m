function p = cat(dim,varargin)
% Concatenate polynomial arrays along specified dimension.
%
% Overwriting matlab.mixin.indexing.RedefinesParen.cat

if nargin < 2
    % no polynomials to concatenate
    p = casos.PS;
    return

elseif length(varargin) == 1
    % single input
    p = varargin{1};
    return
end

% else
assert(any(dim==[1 2]),'Operating dimension must be either 1 or 2.')

% number of polynomials to concatenate
m = length(varargin);

% iterate over polynomials once
% get input parameters, degrees, and variables
[szi,nvi,indets] = cellfun(@(q) getparam(q), varargin, 'UniformOutput',false);

% check dimensions
sizes = vertcat(szi{:});

assert(min(sizes(:,3-dim)) == max(sizes(:,3-dim)), 'Dimensions of polynomials being concatenated are not consistent.')

p = casos.PS;

% size of output
sz(dim) = sum(sizes(:,dim)); sz(3-dim) = sizes(1,3-dim);

% combine variables
[indets,~,ic] = unique([indets{:}]);

% iterate over polynomials twice
% extend degree matrices to combined variables
% extend coefficients to new dimensions
[cfi,dgi] = arrayfun(@(i) extendDC(i,varargin{i},dim,sizes,sz,nvi,ic), 1:m, 'UniformOutput',false);

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(vertcat(cfi{:}), vertcat(dgi{:}));

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,indets);

% new polynomial
p.coeffs = coeffs;
p.degmat = degmat;
p.indets = indets;
p.matdim = sz;

end

function [sz,nv,indeti] = getparam(q)
% Return size, nterms, nvars, nonzero degrees, and indeterminates.

    q = casos.PS(q);

    sz = size(q);
    nv = q.nvars;
%     [i,j,d] = find(q.degmat);
    indeti = q.indets;
end

function [cf,dg] = extendDC(i,q,dim,sizes,sz,nvi,ic)
% Return extended degree and coefficient matrices for i-th polynomial.

    q = casos.PS(q);

    % get nonzero degrees
    [ii,ji,di] = find(q.degmat);

    % extend degree matrix
    dg = sparse(ii,ic(sum([nvi{1:i-1}])+ji),di,q.nterm,max(ic));

    % extend coefficient matrix
    % let a1,...,aL, b1,...,bM, c1,...,cN be the columns of A, B, C
    coeffs = q.coeffs;
    idx = find(sparsity(coeffs));
    [ri,ci] = ind2sub(size(coeffs),idx);
    switch (dim)
        case 1
            % vertical concatenation
            % vec([A;B;C]) = [(a1 b1 c1) ... (aL bM cN)], assuming L=M=N

            % horizontal offset of i-th polynomial's 1st column
            offset = sum(sizes(1:i-1,1)) + (sz(1)-sizes(i,1))*floor((ci-1)/sizes(i,1));

        case 2
            % horizontal concatenation
            % vec([A B C]) = [a1...aL b1...bM c1...cN]

            % horizontal offset of i-th polynomial's coefficient matrix
            offset = sum(sizes(1:i-1,1).*sizes(1:i-1,2));
    end

    % sparsity pattern of extended coefficient matrix (casadi uses 0-index)
    S = casadi.Sparsity.triplet(q.nterm,prod(sz),ri-1,offset+ci-1);
    % assign i-th polynomial's coefficients to extended matrix
    cf = casadi.SX(S,coeffs(idx));
end
