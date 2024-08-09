function b = mpower(a,n)
% Power of polynomial matrix.

assert(size(a,1) == size(a,2), 'Matrix must be square.')
assert(isscalar(n), 'Exponent must be scalar.')

if isscalar(a)
    % fall back to scalar power.
    b = power(a,n);
    return

elseif n == 0
    % identity matrix
    b = a.eye(length(a));
    return

elseif n == 1
    % no change
    b = a;
    return
end

b = a.new_poly;

% Call mpower recursively
[S,b.coeffs] = mpower_internal(a.get_sparsity,a.coeffs,n);

% Note: A for loop on mtimes requires n calls to mtimes. A recursion
% only requires roughly log2(n) calls to mtimes.

b = set_sparsity(b,S);

end
