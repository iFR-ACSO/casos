function tf = is_wellposed(obj)
% Check if polynomial is well posed.

tf = is_wellposed@casos.package.core.GenericPolynomial(obj) ... check sparsity
    && is_equal(sparsity(obj.coeffs), coeff_sparsity(obj.get_sparsity));

end
