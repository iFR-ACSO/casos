function c = plus(a,b)
% Add (join) two polynomial sparsity patterns.

a = casos.Sparsity(a);
b = casos.Sparsity(b);

% prepare error message for incompatible sizes
errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';

% sparsity patterns must be of same size
assert(isequal(size(a),size(b)), errsz, size2str(a), size2str(b))

% join coefficient matrices
c = coeff_plus(a,b,a.coeffs,b.coeffs);

end
