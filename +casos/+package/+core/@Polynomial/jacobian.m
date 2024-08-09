function D = jacobian(f,x)
% Compute symbolic Jacobian matrix of vector polynomial expression.
%
% Use of JACOBIAN for the polynomial Jacobian matrix is not allowed. 
% Use NABLA instead.

if is_indet(x)
    % not allowed: compute polynomial Jacobian matrix.
    error('Use NABLA for the polynomial Jacobian.')
end

% else:
assert(is_symbolic(x),'Second argument must be symbolic polynomial.')

% project x to basis
[X,zi] = poly2basis(x);

% project f to basis
[F,zo] = poly2basis(f);

% G is jacobian of coefficients
G = jacobian(F,X);

% return operator
D = casos.package.operator(G,zi,zo);

end
