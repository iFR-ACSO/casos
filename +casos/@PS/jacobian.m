function g = jacobian(f,x)
% Compute symbolic Jacobian matrix of vector polynomial expression.
%
% Use of JACOBIAN for the polynomial Jacobian matrix is deprecated. 
% Use NABLA instead.

if is_indet(x)
    % DEPRECATED: compute polynomial Jacobian matrix.
    warning('DEPRECATED: Use NABLA for the polynomial Jacobian.')
    g = nabla(f,x);
    return
end

% else
assert(is_symbolic(x),'Second argument must be symbolic polynomial.')

end
