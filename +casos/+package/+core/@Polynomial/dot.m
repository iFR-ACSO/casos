function r = dot(p,q)
% Scalar dot product of two polynomials.

assert(numel(p) == numel(q),'Inputs must be of compatible size.')

% expand coefficient matrices to match
[cf1,cf2] = coeff_expand(p.get_sparsity,q.get_sparsity,p.coeffs,q.coeffs);

% return dot product of coefficients
r = dot(cf1,cf2);

end
