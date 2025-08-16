function [S,coeffs] = coeff_repmat(obj,coeffs,rep)
% Repeat copies of polynomial.

assert(isrow(rep) && length(rep) == 2, 'Replication factors must be a pair (row) of integers or two integer scalars.')

% new dimensions
new_sz = rep.*size(obj);

% repeat coefficients
coeffs = reshape(repmat(coeffs,rep),obj.nterm,prod(new_sz));

% new sparsity pattern
S = new_from_coefficients(coeffs,obj.degmat,obj.indets,new_sz);

end
