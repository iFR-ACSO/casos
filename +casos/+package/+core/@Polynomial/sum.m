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
    B = A.zeros(sizeofMatrixOp(A.get_sparsity,dim));
    return
end

% else
B = A.new_poly;

% finish operation on coefficients
[S,B.coeffs] = coeff_sum(A.get_sparsity,A.coeffs,dim);

% set sparsity
B = set_sparsity(B,S);

end
