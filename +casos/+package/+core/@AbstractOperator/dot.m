function c = dot(a,b)
% Dot operation on linear operators.
%
% If b is another operator, then dot(a,b) is a composition.
% If b is a polynomial, then dot(a,b) is an evaluation.

if isa(b,'casos.package.core.GenericPolynomial')
    % evaluation
    assert(isequal(size_in(a), size(b)),'Input dimensions mismatch.')

    B = poly2basis(b,a.sparsity_in);
    % evaluate operator on nonzero coordinates
    C = a.matrix*B;
    % return polynomial
    c = casos.package.polynomial(a.sparsity_out,C);

elseif isa(b,'casos.package.core.AbstractOperator')
    % composition
    assert(isequal(size_in(a), size_out(b)),'Dimensions mismatch for composition.')

    % find common nonzero coordinates
    [~,I1,I2] = op_intersect(a.sparsity_in,b.sparsity_out);

    % composite operators
    C = a.matrix(:,I1)*b.matrix(I2,:);
    % return new operator
    c = a.new_operator(C,a.sparsity_out,b.sparsity_in);

else
    error('Function "dot" not defined for input %s.',class(b))
end

end
