function c = ldivide(a,p)
% Element-wise left array division with constant matrix.

assert(is_zerodegree(a),'Only division by constant or symbolic matrix possible.')

if is_operator(p)
    % divide operator by scalar
    assert(isscalar(a),'Only division by scalar allowed for operators.')

    % divide coefficients by scalar
    coeffs = a.coeffs .\ p.coeffs;
    Sp = p.get_sparsity;

else

% dimensions are compatible if equal or one summand is row/column
if ~check_sz_comptbl(a,p)
    throw(casos.package.core.IncompatibleSizesError.basic(a,p));
end

% input dimensions
sza = size(a);
szp = size(p);

% dimensions of element-wise product
sz = max(sza,szp);

% find zero dimension
I0 = (sza == 0) | (szp == 0);

% handle simple case(s) for speed up
if isempty(a) || isempty(p)
    % element-wise multiplication with empty polynomial is empty
    sz(I0) = 0;

    c = p.empty(sz);
    return
end

% else
c = p.new_poly;

% reshape to output dimensions
[Sa,cfa] = coeff_repmat(a.get_sparsity,a.coeffs,sz./sza);
[Sp,cfp] = coeff_repmat(p.get_sparsity,p.coeffs,sz./szp);

% element-wise division of coefficients
coeffs = coeff_repterms(Sa,cfa,p.nterm) .\ cfp;

end

% update sparsity pattern
[S,c.coeffs] = coeff_update(Sp,coeffs);

% set sparsity
c = set_sparsity(c,S);

end
