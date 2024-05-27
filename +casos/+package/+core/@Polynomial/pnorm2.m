function r = pnorm2(a)
% Compute polynomial 2-norm as integral of a^2 over [0 1].
%
% Note: This 2-norm does NOT correspond to the norm induced by dot.

% transpose
[St,cft] = coeff_transpose(a.get_sparsity,a.coeffs);

% TODO: use internal operations
[Sb,cfb] = mtimes_internal(St,a.get_sparsity,cft,a.coeffs);

% integral over [0 1]
r = coeff_properint(Sb,cfb);

end
