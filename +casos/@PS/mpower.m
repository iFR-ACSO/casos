function b = mpower(a,n)
% Power of polynomial matrix.

assert(size(a,1) == size(a,2), 'Matrix must be square.')
assert(all(size(n) <= 1), 'Matrix-valued exponent not supported.')
assert(length(a) == 1 || ~isempty(n), 'Incorrect dimensions.')

if length(a) <= 1
    % fall back to scalar power.
    b = power(a,n);
    return

elseif n == 0
    % identity matrix
    b = casos.PS.eye(length(a));
    return

elseif n == 1
    % no change
    b = a;
    return
end

% Call mpower recursively
% from multipoly
% Note: A for loop on mtimes requires n calls to mtimes. A recursion
% only requires roughly log2(n) calls to mtimes.
if floor(n/2)==ceil(n/2)
    b = mpower(a,n/2);
    b = mtimes(b,b);

else
    b = mpower(a,(n-1)/2);
    b = mtimes(b,b);
    b = mtimes(b,a);
end
% TODO: don't remove zero coefficients and/or degrees within recursion?

end