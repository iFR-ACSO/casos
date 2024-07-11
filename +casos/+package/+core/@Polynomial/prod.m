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
    B = A.ones(sizeofMatrixOp(A.get_sparsity,dim));
    return

elseif isequal(dim,'all') && isscalar(A) ...
        || ~isequal(dim,'all') && size(A,dim) == 1
    % nothing to do
    B = A;
    return
end

% else:
B = A.new_poly;

% compute product of coefficients
[S,B.coeffs] = coeff_prod(A.get_sparsity,A.coeffs,dim);

% set sparsity
B = set_sparsity(B,S);

end
