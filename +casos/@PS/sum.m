function B = sum(A,dim)
%SUM Overloads Matlabs internal sum function to work with polynomials of type casos.PS
% Inputs:
%   A: Polynomial matrix of arbitrary size
%   dim (default='all'): Defines wether columns, rows or all elements are summed up
%       dim = 1: sum along first dimension (rows)
%       dim = 2: sum along second dimension (columns)
%       dim = 'all': sum all elements
% Outputs:
%   B: matrix of polynomials
%
% Example with matrix of polynomials:
% p1, p2, p3, p4 are casos.PS polynomials (with symbolic or numeric coefficients)
% A = [p1 p2; p3 p4];
% sum(A,'all')      outputs i.e.                    [p1+p2+p3+p4]
% sum(A,1)          outputs sum along columns i.e.  [p1+p3, p2+p4]
% sum(A,2)          outputs sum along row i.e.      [p1+p2; p3+p4]

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
evalMatrixOp(cfa,@sum,dim);

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
