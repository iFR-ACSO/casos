function c = mtimes(a,b)
% Matrix multiplication of two polynomials.

assert(~is_operator(a) && ~is_operator(b), 'Not allowed for operators.')

if isempty(a) || isempty(b)
    % empty multiplication
    c = a.zeros(size(a,1),size(b,2));
    return

elseif isscalar(a) || isscalar(b)
    % fall back to scalar multiplication
    c = times(a,b);
    return

elseif ~check_sz_mtimes(a,b)
    % dimensions are compatible if size(a,2) == size(b,1)
    throw(casos.package.core.IncompatibleSizesError.matrix(a,b));
end

% else
c = a.new_poly;

% compute coefficient matrix and sparsity
[S,c.coeffs] = mtimes_internal(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
