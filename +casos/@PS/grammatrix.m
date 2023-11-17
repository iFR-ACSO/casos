function [Q,Z,K,z] = grammatrix(p,I)
% Attempts to compute a symbolic Gram matrix for polynomial vector.

if nargin < 2
    I = true(size(p));
end
if ~any(I)
    % empty polynomial
    Q = casadi.SX;
    [Z,K,z] = grambasis(casos.PS);

    return
end

% compute Gram basis vector
[Z,K,z] = grambasis(p,I);

% get symbolic variables from coefficients
syms = symvar(p.coeffs(:,find(I)));

% check if dimensions are compatible
if length(syms) == sum(K.^2)

% stack symbols horizontally
Q = horzcat(syms{:});

% attempt to build Gram form
pgram = casos.PS.zeros(size(p));
pgram(I) = (Q*Z)';

% check if Gram form is equal
diff = (p - pgram);
if is_zero(diff.coeffs(:,find(I)))
    % return Gram matrix
    return
end

end

% else:
% return empty Gram matrix
Q = [];

end
