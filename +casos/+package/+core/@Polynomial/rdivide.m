function c = rdivide(p,b)
% Element-wise right array division with constant matrix.

assert(is_zerodegree(b),'Only division by constant or symbolic matrix possible.')

% input dimensions
szp = size(p);
szb = size(b);

% find zero dimension
I0 = (szp == 0) | (szb == 0);

% dimensions are compatible if equal or one summand is row/column
if ~check_sz_comptbl(p,b)
    throw(casos.package.core.IncompatibleSizesError.basic(p,b));
end

% dimensions of element-wise product
sz = max(szp,szb);

% handle simple case(s) for speed up
if isempty(p) || isempty(b)
    % element-wise multiplication with empty polynomial is empty
    sz(I0) = 0;

    c = p.empty(sz);
    return
end

% else
c = p.new_poly;

% reshape to output dimensions
[Sp,cfp] = coeff_repmat(p.get_sparsity,p.coeffs,sz./szp);
[Sb,cfb] = coeff_repmat(b.get_sparsity,b.coeffs,sz./szb);

% element-wise division of coefficients
coeffs = cfp ./ coeff_repterms(Sb,cfb,p.nterm);

% update sparsity pattern
[S,c.coeffs] = coeff_update(Sp,coeffs,sz);

% set sparsity
c = set_sparsity(c,S);

end
