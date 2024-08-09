function H = hessian(f,x)
% Compute symbolic Hessian matrix of scalar, zero-degree polynomial
% expression.
%
% Polynomial Hessian matrix is not supported.

if is_indet(x)
    % not allowed
    warning('Polynomial Hessian matrix is not supported.')
end

assert(is_symbolic(x),'Second argument must be symbolic polynomial.')
assert(isscalar(f) && is_zerodegree(f),'First argument must be scalar, zero-degree polynomial.')

% project x to basis
[X,zi] = poly2basis(x);

% convert f to SX
F = poly2basis(f);

% hessian of coefficients
H = casos.package.operator(hessian(F,X),zi,zi);

end
