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

elseif size(A,dim) == 1
    % nothing to do
    B = A;
    return
end

% else:
B = casos.PS;

nt = A.nterm;
nv = A.nvars;

% reshape input coefficient matrix
[cfa,sz] = prepareMatrixOp(A.coeffs,A.matdim,dim);

% multiply entries
L = numel(A)/prod(sz);
nb = prod(sz);
% (sum_a c_a[1]*x^a).* ... .* (sum_a c_a[L]*x^a)
%  = sum_a1 ... sum_aL (c_a1[1].* ... .* c_aL[L])*x^(a1+...+aL)
I = cell(1,L);
[I{:}] = ndgrid(1:nt); I = cellfun(@(c) c(:),I,'UniformOutput',false);
idx = horzcat(I{:});

% permute degree matrix
%
%       | dax ... day |      | da1x ... da1y        |
%   D = |  :       :  |  ->  |   :        :     ... |
%       | dzx ... dzy |      | daLx ... daLy        |
%
dga = reshape(A.degmat(idx',:),L,nv*nt^L);
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
        ind = sub2ind(size(cfa),row,col);
        C = reshape(cfa(ind),nt^L,L);

        % combine coefficients
        cfb = sx_prod(C,2);

    case 2
        % compute product along second dimension
        % rows are rows of c_a, sorted by term
        row = kron((0:nb-1)',nt*ones(nt^L,L)) + repmat(idx,nb,1);
        col = repmat(1:L,nb*nt^L,1);    % select row
        ind = sub2ind(size(cfa),row,col);
        C = reshape(cfa(ind),nb*nt^L,L);

        % combine coefficients
        cfb = sx_prod(C,2);

    case 1
        % compute product along first dimension
        % columns are colums of c_a, sorted by coefficient
        term0 = nb*(0:nt-1)';
        row = repmat(1:L,nb*nt^L,1);    % select column
        col = kron(term0(idx),ones(nb,1)) + repmat((1:nb)',nt^L,L);
        ind = sub2ind(size(cfa),row,col);
        C = reshape(cfa(ind'),L,nb*nt^L);

        % combine coefficients
        cfb = sx_prod(C,1);
end

% reshape output coefficient matrix
coeffs = finishMatrixOp(cfb,sz,dim);
% reshape output degree matrix
degmat = reshape(dgb,nt^L,nv);

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs, degmat);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,A.indets);

% new polynomial
B.coeffs = coeffs;
B.degmat = degmat;
B.indets = indets;
B.matdim = sz;

end
