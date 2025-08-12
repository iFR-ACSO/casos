function c = plus(a,b)
% Add two polynomials.

assert(is_operator(a) == is_operator(b), 'Must not mix polynomials and operators.')

% input dimensions
sza = size(a);
szb = size(b);

% find zero dimension
I0 = (sza == 0) | (szb == 0);

% dimensions are compatible if equal or one summand is row/column
if ~check_sz_comptbl(a,b)
    throw(casos.package.core.IncompatibleSizesError.basic(a,b));
end

% dimensions of sum
sz = max(sza,szb);

% handle simple case(s) for speed up
if isempty(a) || isempty(b)
    % addition with empty polynomial is empty
    sz(I0) = 0;

    c = a.empty(sz);
    return
end

% else
c = a.new_poly;

% reshape to output dimensions
[S1,cfa] = coeff_repmat(a.get_sparsity,a.coeffs,sz./sza);
[S2,cfb] = coeff_repmat(b.get_sparsity,b.coeffs,sz./szb);

% add coefficient matrices
[S,c.coeffs] = coeff_plus(S1,S2,cfa,cfb);

% set sparsity
c = set_sparsity(c,S);

end
