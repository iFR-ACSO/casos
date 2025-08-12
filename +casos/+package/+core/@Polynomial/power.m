function b = power(a,n)
% Element-wise powers.

assert(~is_operator(a), 'Not allowed for operators.')

% input dimensions
sza = size(a);
szn = size(n);

% find zero dimension
I0 = (sza == 0) | (szn == 0);

% dimensions are compatible if equal or one summand is row/column
if ~check_sz_comptbl(a,n)
    throw(casos.package.core.IncompatibleSizesError.basic(a,n));
end

% dimensions of sum
sz = max(sza,szn);

% handle simple case(s) for speed up
if isempty(a) || isempty(n)
    % element-wise power with empty polynomial/exponent is empty
    sz(I0) = 0;

    b = a.new_poly(sz);
    return

elseif all(n==0)
    % power of zero is one
    b = a.ones(sz);
    return
end
% TODO: handle or escape for other simple cases, e.g., scalar, constant
% matrix, single term etc.?

% else
b = a.new_poly;

% reshape to output dimensions
[S,coeffs] = coeff_repmat(a.get_sparsity,a.coeffs,sz./sza);
deg = repmat(n,sz./szn);

% compute power of coefficient matrix
[S,b.coeffs] = coeff_power(S,coeffs,deg);

% set sparsity
b = set_sparsity(b,S);

end
