function [G,zi,zo] = jacobian(f,x)
% Compute symbolic Jacobian matrix of vector polynomial expression.
%
% Use of JACOBIAN for the polynomial Jacobian matrix is deprecated. 
% Use NABLA instead.

if is_indet(x)
    % DEPRECATED: compute polynomial Jacobian matrix.
    warning('DEPRECATED: Use NABLA for the polynomial Jacobian.')
    G = nabla(f,x);
    return
end

% else
f = casos.PS(f);
x = casos.PS(x);

assert(is_symbolic(x),'Second argument must be symbolic polynomial.')

% project to basis
[X,zi] = poly2basis(x);
[F,zo] = poly2basis(f);

% G is jacobian of coefficients
G = jacobian(F,X);

end
