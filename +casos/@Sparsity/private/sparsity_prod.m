function coeffs = sparsity_prod(coeffs,sz,dim)
% Compute product of coefficient sparsity pattern along given dimension.

cfa = pattern_inverse(coeffs);

% compute sum of inverse
cfb = sparsity_sum(cfa,sz,dim);

% product is inverse of sum
coeffs = pattern_inverse(cfb);

end
