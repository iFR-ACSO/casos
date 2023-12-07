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

% sort to match order in basis
[ii,jj] = ind2sub(size(Z.coeffs),find(sparsity(Z.coeffs)));
% only sort for each entry in p
nZ = size(Z,1);
[~,idx] = sort(ii + ceil(jj/nZ)*nZ);

% stack symbols horizontally
Q = casadi.SX(length(syms),1);
Q(idx) = horzcat(syms{:});

% attempt to build Gram form
pgram = casos.PS.zeros(size(p));
pgram(I) = casos.PS(Z,Q);

% check if Gram form is equal
diff = (p - pgram);
dcfs = diff.coeffs(:,find(I));
if is_zero(dcfs) || is_zero(simplify(dcfs))
    % return Gram matrix
    return
end

end

% else:
% return empty Gram matrix
Q = [];

end
