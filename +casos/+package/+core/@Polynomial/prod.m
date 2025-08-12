function b = prod(a,dim)
% Return product of array elements.

assert(~is_operator(a), 'Not allowed for operators.')

if nargin > 1
    % nothing to do
elseif isempty(a) && ~isvector(a)
    % product of empty matrix (return scalar)
    dim = 'all';
elseif isrow(a)
    % product along second dimension (return column)
    dim = 2;
else
    % product along first dimension (return row)
    dim = 1;
end

if isempty(a)
    % product of empty matrix is one
    b = a.ones(sizeofMatrixOp(a.get_sparsity,dim));
    return

elseif isequal(dim,'all') && isscalar(a) ...
        || ~isequal(dim,'all') && size(a,dim) == 1
    % nothing to do
    b = a;
    return
end

% else:
b = a.new_poly;

% compute product of coefficients
[S,b.coeffs] = coeff_prod(a.get_sparsity,a.coeffs,dim);

% set sparsity
b = set_sparsity(b,S);

end
