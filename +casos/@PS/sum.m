function B = sum(A,dim)
% Return sum of array elements.

if nargin > 1
    % nothing to do
elseif isempty(A) && ~isvector(A)
    % sum of empty matrix (return scalar)
    dim = 'all';
elseif isrow(A)
    % sum along second dimension (return column)
    dim = 2;
else
    % sum along first dimension (return row)
    dim = 1;
end

if isempty(A)
    % sum of empty matrix is zero
    [~,sz] = prepareMatrixOp([],A.matdim,dim);

    B = casos.PS.zeros(sz);
    return
end

% else
B = casos.PS;

% reshape input coefficient matrix
[cfa,sz] = prepareMatrixOp(A.coeffs,A.matdim,dim);

% sum coefficients
switch (dim)
    case 'all'
        % compute total sum of each coefficient
        cfb = sum(cfa,2);
    otherwise
        % compute sum of coefficient along dimension
        cfb = sum(cfa,dim);
end

% reshape output coefficient matrix
coeffs = finishMatrixOp(cfb,sz,dim);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,A.degmat,A.indets);

% new polynomial
B.coeffs = coeffs;
B.degmat = degmat;
B.indets = indets;
B.matdim = sz;

end
