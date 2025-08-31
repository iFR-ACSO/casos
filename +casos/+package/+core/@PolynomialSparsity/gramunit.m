function Z = gramunit(w)
% Return a Gram unit basis for a given monomial sparsity pattern.

% double the exponents of each monomial
Z = new_from_coefficients(w.coeffs,2*w.degmat,w.indets,w.matdim);

end
