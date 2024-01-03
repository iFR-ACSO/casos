function c = ldivide(a,p)
% Element-wise left array division with constant or symbolic matrix.

p = casos.PS(p);
a = casos.PS(a);

assert(is_zerodegree(a),'Only division by constant or symbolic matrix possible.')

% input dimensions
sza = size(a);
szp = size(p);

% prepare error message for incompatible sizes
errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';

% compare dimensions
I = (szp == sza);
% find zero dimension
I0 = (szp == 0) | (sza == 0);
% find one dimension
I1 = (szp == 1) | (sza == 1);

% dimensions are compatible if equal or one factor is row/column
assert(all(I | I1), errsz, size2str(sza), size2str(szp))

% dimensions of element-wise product
sz = max(szp,sza);

% handle simple case(s) for speed up
if isempty(p) || isempty(a)
    % element-wise multiplication with empty polynomial is empty
    sz(I0) = 0;

    c = casos.PS.zeros(sz);
    return
end

% else
c = casos.PS;

nt = p.nterm;

% reshape to output dimensions
cfa = reshape(repmat(a.coeffs,sz./sza), 1,prod(sz));
cfp = reshape(repmat(p.coeffs,sz./szp),nt,prod(sz));

% a.\(sum_b c_b*x^b) = sum_a (a.\c_b)*x^b
coeffs = repmat(cfa,nt,1) .\ cfp;

% new polynomial
c.coeffs = coeffs;
c.degmat = p.degmat;
c.indets = p.indets;
c.matdim = sz;

end
