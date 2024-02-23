function Z = adjoint_inverse(obj)
% Compute the adjoint-inverse of a Gram basis.

scale = sum(obj.coeffs,2);

% new polynomial
Z = casos.PS;
Z.coeffs = scale.\obj.coeffs;
Z.degmat = obj.degmat;
Z.indets = obj.indets;
Z.matdim = obj.matdim;

end
