function tf = is_linear(p,q)
% Check if polynomial is linear in symbols.

p = casos.PS(p);
q = casos.PS(q);

% check if polynomial is in symbolic Gram form
Q = grammatrix(q);

if ~isempty(Q)
    tf = is_linear(p.coeffs,Q(:));
    return
end

% else
assert(is_symbolic(q),'Second argument must be purely symbolic.')

tf = is_linear(p.coeffs,q.coeffs);

end
