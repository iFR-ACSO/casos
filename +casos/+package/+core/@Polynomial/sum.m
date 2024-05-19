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

% prepare for operation on coefficients
[S,coeffs] = prepareMatrixOp(A.get_sparsity,A.coeffs,dim);

if isempty(A)
    % sum of empty matrix is zero
    B = A.new_poly(size(S));
    return
end

% else
B = A.new_poly;

% sum coefficients
cfb = evalMatrixOp(coeffs,@sum,dim);

% finish operation on coefficients
[S,B.coeffs] = coeff_sum(S,cfb,dim);

% set sparsity
B = set_sparsity(B,S);

end
