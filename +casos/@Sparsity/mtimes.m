function c = mtimes(a,b)
% Multiply (intersect) two polynomial sparsity patterns.

a = casos.Sparsity(a);
b = casos.Sparsity(b);

% sparsity patterns must be of same size
if ~check_sz_equal(a,b)
    throw(casos.package.core.IncompatibleSizesError.other(a,b));
end

% intersect coefficient matrices
c = op_intersect(a,b);

end
