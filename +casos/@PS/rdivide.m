function c = rdivide(p,b)
% Element-wise right array division with constant or symbolic matrix.

p = casos.PS(p);
b = casos.PS(b);

assert(is_zerodegree(b),'Only division by constant or symbolic matrix possible.')

% input dimensions
szp = size(p);
szb = size(b);

% prepare error message for incompatible sizes
errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';

% compare dimensions
I = (szp == szb);
% find zero dimension
I0 = (szp == 0) | (szb == 0);
% find one dimension
I1 = (szp == 1) | (szb == 1);

% dimensions are compatible if equal or one factor is row/column
assert(all(I | I1), errsz, size2str(szp), size2str(szb))

% dimensions of element-wise product
sz = max(szp,szb);

% handle simple case(s) for speed up
if isempty(p) || isempty(b)
    % element-wise multiplication with empty polynomial is empty
    sz(I0) = 0;

    c = casos.PS.zeros(sz);
    return
end

% else
c = casos.PS;

nt = p.nterm;

% reshape to output dimensions
cfp = reshape(repmat(p.coeffs,sz./szp),nt,prod(sz));
cfb = reshape(repmat(b.coeffs,sz./szb), 1,prod(sz));

% (sum_a c_a*x^a)./b = sum_a (c_a./b)*x^a
coeffs = cfp ./ repmat(cfb,nt,1);

% new polynomial
c.coeffs = coeffs;
c.degmat = p.degmat;
c.indets = p.indets;
c.matdim = sz;

end
