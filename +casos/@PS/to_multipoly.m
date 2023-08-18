function q = to_multipoly(p)
% Convert a polynomial with constant coefficients into a multipoly object.

assert(~is_symexpr(p), 'Cannot convert polynomials with symbolic coefficients.')

% get coefficients
coeffs = casadi.DM(p.coeffs);

q = polynomial(sparse(coeffs),p.degmat,p.indets,p.matdim);

end
