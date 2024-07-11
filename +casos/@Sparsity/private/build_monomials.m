function z = build_monomials(degmat,indets)
% Build a monomial sparsity pattern.

z = casos.Sparsity;

% number of monomials
nt = size(degmat,1);

% set output
z.coeffs = casadi.Sparsity.dense(nt,1);
z.degmat = sparse(degmat);
z.indets = indets;
z.matdim = [1 1];

end
