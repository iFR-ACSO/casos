function r = dot(p,q)
% Scalar dot product of two polynomials.

p = casos.PS(p);
q = casos.PS(q);

assert(numel(p) == numel(q),'Inputs must be of compatible size (LHS has %g elements, RHS has %g).',numel(p),numel(q))

% get basis of first argument
[P,Z] = poly2basis(p);
% project second argument
Q = poly2basis(q,Z);

% return dot product of coefficients
r = dot(P,Q);

end
