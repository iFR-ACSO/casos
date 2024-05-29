function M = expand(matrix,sz,Iin,Iout)
% Expand operator matrix to joint nonzeros.

[ii,jj] = get_triplet(sparsity(matrix));

% new sparsity pattern of operator matrix
SpA = casadi.Sparsity.triplet(sz(1),sz(2),Iout(ii+1)-1,Iin(jj+1)-1);

% return expanded operator matrix
M = sparsity_cast(matrix,SpA);

end
