function c = plus(a,b)
% Add two polynomials.

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

% dimensions are compatible if equal or one summand is row/column
assert(all(I | I1), errsz, size2str(a), size2str(b))

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