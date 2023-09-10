function [Q,Z,K,z] = grammatrix(p)
% Attempts to compute a symbolic Gram matrix for polynomial vector.

% compute Gram basis vector
[Z,K,z] = grambasis(p);

% get symbolic variables from coefficients
syms = symvar(p.coeffs);

% check dimensions if are compatible
if length(syms) == sum(K.^2)

% stack symbols horizontally
Q = horzcat(syms{:});

% attempt to build Gram form
pgram = (Q*Z)';

% check if Gram form is equal
if isequal(p, pgram)
    % return Gram matrix
    return
end

end

% else:
% return empty Gram matrix
Q = [];

end
