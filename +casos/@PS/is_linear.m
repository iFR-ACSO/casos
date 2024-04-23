function tf = is_linear(p,q,I)
% Check if polynomial is linear in symbols.

p = casos.PS(p);
q = casos.PS(q);

if nargin < 3
    I = true(size(p));
end

% check if polynomial is in symbolic Gram form
Q = grammatrix(q);

if ~isempty(Q)
    tf = is_linear(p.coeffs(:,find(I)),Q(:));
    return
end

% else
[is_sym,coeffs] = is_symbolic(q);
assert(is_sym,'Second argument must be purely symbolic.')

tf = is_linear(p.coeffs(:,find(I)),coeffs);

end
