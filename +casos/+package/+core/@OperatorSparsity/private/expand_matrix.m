function M = expand_matrix(S_coeffs,coeffs,sz,Iin,Iout)
% Expand operator matrix to joint nonzeros.

[ii,jj] = get_triplet(S_coeffs);

% new sparsity pattern of operator matrix
SpA = casadi.Sparsity.triplet(sz(1),sz(2),Iout(ii+1)-1,Iin(jj+1)-1);

% return expanded operator matrix
M = sparsity_cast(coeffs,SpA);

end
