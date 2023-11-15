function z = build_monomials(degmat,indets)
% Build a monomial vector.

z = casos.PS;

% number of monomials
nt = size(degmat,1);

% set output
z.coeffs = casadi.SX.eye(nt);
z.degmat = sparse(degmat);
z.indets = indets;
z.matdim = [nt 1];

end
