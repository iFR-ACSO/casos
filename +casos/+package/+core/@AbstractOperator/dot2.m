function r = dot2(obj,q,p)
% Evaluate operator as bilinear product.
%
% This is equivalent to dot(q,op(p)).

assert(isequal(size_in(obj), size(p)),'Input dimensions mismatch.')
assert(isequal(size_out(obj), size(q)),'Input dimensions mismatch.')

% project onto lhs and rhs sparsity patterns
P = poly2basis(p,obj.sparsity_in);
Q = poly2basis(q,obj.sparsity_out);

% evaluate operator on nonzero coordinates
r = Q'*obj.matrix*P;

end
