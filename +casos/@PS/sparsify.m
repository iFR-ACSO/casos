function b = sparsify(a)
% Sparsify symbolic expressions.

coeffs = sparsify(a.coeffs);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,a.degmat,a.indets);

% new polynomial
b = casos.PS;
b.coeffs = coeffs;
b.degmat = degmat;
b.indets = indets;
b.matdim = a.matdim;

end
