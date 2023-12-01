function tf = is_linear(obj,s)
% Check if polynomial is linear in symbols.

s = casos.PS(s);

% check if polynomial is in symbolic Gram form
Q = grammatrix(s);

if ~isempty(Q)
    tf = is_linear(obj.coeffs,Q(:));
    return
end

% else
assert(is_symbolic(s),'Second argument must be purely symbolic.')

tf = is_linear(obj.coeffs,s.coeffs);

end
