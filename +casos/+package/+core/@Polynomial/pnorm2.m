function r = pnorm2(a)
% Compute polynomial 2-norm as integral of a^2 over [0 1].
%
% Note: This 2-norm does NOT correspond to the norm induced by dot.

% reshape to vector
a = reshape(a,numel(a),1);
Sa = a.get_sparsity;

% TODO: use internal operations
[Sb,cfb] = mtimes_internal(T(Sa),Sa,a.coeffs,a.coeffs);

% integral over [0 1]
r = coeff_properint(Sb,cfb);

end
