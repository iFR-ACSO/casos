function q = to_multipoly(p)
% Convert a polynomial with constant coefficients into a multipoly object.

assert(is_constant(p), 'Cannot convert polynomials with symbolic coefficients.')

p = casos.PD(p);

% get coefficients
c = list_of_coeffs(p);
d = list_of_degree(p);
v = list_of_indets(p);

coeffs = reshape(horzcat(c{:}),numel(p),p.nterm)';
degmat = vertcat(d{:});

q = polynomial(sparse(coeffs),degmat,v,size(p));

end
