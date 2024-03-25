function [G,zi,zo] = jacobian(f,x)
% Compute symbolic Jacobian matrix of vector polynomial expression.
%
% Use of JACOBIAN for the polynomial Jacobian matrix is deprecated. 
% Use NABLA instead.

x = casos.PS(x);

if is_indet(x)
    % DEPRECATED: compute polynomial Jacobian matrix.
    warning('DEPRECATED: Use NABLA for the polynomial Jacobian.')
    G = nabla(f,x);
    return
end

% else:
% check x for Gram symbolic
[X,zi] = grammatrix(x);

if isempty(X)
    % polynomial is not Gram
    assert(is_symbolic(x),'Second argument must be symbolic polynomial.')

    % project x to basis
    [X,zi] = poly2basis(x);
end

% project f to basis
f = casos.PS(f);
[F,zo] = poly2basis(f);

% G is jacobian of coefficients
G = jacobian(F,X);

end
