function tf = is_linear(p,q)
% Check if polynomial is linear in symbols.

assert(is_symbolic(q),'Second argument must be purely symbolic.')

% get nonzero coordinates
Q = coordinates(q);

% check if coefficients are linear
tf = is_linear(p.coeffs,Q);

end
