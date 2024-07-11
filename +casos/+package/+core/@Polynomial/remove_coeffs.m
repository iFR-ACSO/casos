function obj = remove_coeffs(obj,tol)
% Remove coefficients that are below a given tolerance.

% identify coefficients to keep
S_coeffs = sparsity(sparsify(abs(obj.coeffs) >= tol));

% project onto coefficient sparsity
coeffs = project(obj.coeffs,S_coeffs);

% update coefficients
[S,obj.coeffs] = coeff_update(obj.get_sparsity,coeffs);

% set sparsity
obj = set_sparsity(obj,S);

end
