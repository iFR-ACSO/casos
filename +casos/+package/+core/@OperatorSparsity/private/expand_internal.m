function [cf1,cf2,Si,So] = expand_internal(obj,S2,cfa,cfb)
% Internal function for expansion of operator coefficient matrices.

% join input nonzero coordinates
[Si,I1i,I2i] = op_join(obj.sparsity_in,S2.sparsity_in);
% join output nonzero coordinates
[So,I1o,I2o] = op_join(obj.sparsity_out,S2.sparsity_out);

% size of resulting operator matrix
sz = [nnz(So) nnz(Si)];

% expand matrices to joint nonzeros
cf1 = expand(cfa,sz,I1i,I1o);
cf2 = expand(cfb,sz,I2i,I2o);

end

function M = expand(matrix,sz,Iin,Iout)
% Expand operator matrix to joint nonzeros.

[ii,jj] = get_triplet(sparsity(matrix));

% new sparsity pattern of operator matrix
SpA = casadi.Sparsity.triplet(sz(1),sz(2),Iout(ii+1)-1,Iin(jj+1)-1);

% return expanded operator matrix
M = sparsity_cast(matrix,SpA);

end
