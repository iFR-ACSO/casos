function c = mrdivide(p,a)
% Matrix right division with constant or symbolic matrix.
%
% Note: c = p/A is equivalent to solving the linear equation c*A = p.

p = casos.PS(p);
a = casos.PS(a);

assert(is_zerodegree(a),'Only division by constant or symbolic matrix possible.')

if isscalar(a)
    % fall back to scalar division
    c = rdivide(p,a);
    return
end

% input dimensions
szp = size(p);
sza = size(a);

% prepare error message for incompatible sizes
errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';

% dimensions are compatible if numbers of columns agree
assert(szp(2) == sza(2), errsz, size2str(szp), size2str(sza))

% else
c = casos.PS;

% reshape coefficients into list of row vectors
cfp = prepareMatrixOp(p.coeffs,szp,2);

% (sum_b c_b*x^b)/A = sum_b (c_a/A)*x^b
cfc = cfp / casadi.SX(a);

% output dimension (see MATLAB mrdivide)
sz = [szp(1) sza(1)];

% reshape output coefficients matrix
coeffs = finishMatrixOp(cfc,sz,2);

% new polynomial
c.coeffs = coeffs;
c.degmat = p.degmat;
c.indets = p.indets;
c.matdim = sz;

end
