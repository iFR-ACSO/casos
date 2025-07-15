function Z = gramunit(w)
% Return a Gram unit basis for a given monomial sparsity pattern.

Z = casos.Sparsity(w);
% double the exponents of each monomial
Z.degmat = 2*w.degmat;

end
