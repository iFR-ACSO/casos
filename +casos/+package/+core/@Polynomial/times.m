function c = times(a,b)
% Element-wise multiplication of two polynomial matrices.

% input dimensions
sza = size(a);
szb = size(b);

% find zero dimension
I0 = (sza == 0) | (szb == 0);

% dimensions are compatible if equal or one summand is row/column
if ~check_sz_comptbl(a,b)
    throw(casos.package.core.IncompatibleSizesError.basic(a,b));
end

% dimensions of element-wise product
sz = max(sza,szb);

% handle simple case(s) for speed up
if isempty(a) || isempty(b)
    % element-wise multiplication with empty polynomial is empty
    sz(I0) = 0;

    c = a.empty(sz);
    return
end
% TODO: handle or escape for other simple cases, e.g., scalar, constant
% matrix, single term etc.?

% else
c = a.new_poly;

% reshape to output dimensions
[S1,cfa] = coeff_repmat(a.get_sparsity,a.coeffs,sz./sza);
[S2,cfb] = coeff_repmat(b.get_sparsity,b.coeffs,sz./szb);

% multiply coefficient matrices
[S,c.coeffs] = coeff_times(S1,S2,cfa,cfb);

% set sparsity
c = set_sparsity(c,S);

end
