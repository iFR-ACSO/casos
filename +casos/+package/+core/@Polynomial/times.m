function c = times(a,b)
% Element-wise multiplication of two polynomial matrices.

% input dimensions
sza = size(a);
szb = size(b);

% prepare error message for incompatible sizes
errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';

% compare dimensions
I = (sza == szb);
% find zero dimension
I0 = (sza == 0) | (szb == 0);
% find one dimension
I1 = (sza == 1) | (szb == 1);

% dimensions are compatible if equal or one factor is row/column
assert(all(I | I1), errsz, size2str(a), size2str(b))

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
