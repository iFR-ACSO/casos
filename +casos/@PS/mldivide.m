function c = mldivide(a,p)
% Matrix left division with constant or symbolic matrix.
%
% Note: c = A\p is equivalent to solving the linear equation A*c = p.

a = casos.PS(a);
p = casos.PS(p);

assert(is_zerodegree(a),'Only division by constant or symbolic matrix possible.')

if isscalar(a)
    % fall back to scalar division
    c = ldivide(a,p);
    return
end

% input dimensions
sza = size(a);
szp = size(p);

% prepare error message for incompatible sizes
errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';

% dimensions are compatible if numbers of rows agree
assert(szp(1) == sza(1), errsz, size2str(sza), size2str(szp))

% else
c = casos.PS;

% reshape coefficients into list of row vectors
cfp = prepareMatrixOp(p.coeffs,szp,1);

% A\(sum_b c_b*x^b) = sum_b (A\c_a)*x^b
cfc = casadi.SX(a) \ cfp;

% output dimension (see MATLAB mrdivide)
sz = [sza(2) szp(2)];

% reshape output coefficients matrix
coeffs = finishMatrixOp(cfc,sz,1);

% new polynomial
c.coeffs = coeffs;
c.degmat = p.degmat;
c.indets = p.indets;
c.matdim = sz;

end
