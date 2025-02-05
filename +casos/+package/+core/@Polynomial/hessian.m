function D = hessian(f,x)
% Compute symbolic Hessian matrix of vector polynomial expression.


if is_indet(x)
    % not allowed: compute polynomial Hessian matrix.
    error('Use NABLA for the polynomial Jacobian.')
end

% else:
assert(is_symbolic(x),'Second argument must be symbolic polynomial.')

% project x to basis
[X,zi] = poly2basis(x);

% project f to basis
[F,zo] = poly2basis(f);

% G is jacobian of coefficients
G = hessian(F,X);

% return operator
D = casos.package.operator(G,zi,zi);

end
