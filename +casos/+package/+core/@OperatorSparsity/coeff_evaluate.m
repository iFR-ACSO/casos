function [S,coeffs] = coeff_evaluate(obj,S2,coeff1,coeff2)
% Compute coefficient matrix for operator evaluation.

assert(numel(obj.sparsity_in) == numel(S2), 'Inputs must be of compatible size.')

% find nonzero input coordinates
[~,I1,I2] = op_intersect(obj.sparsity_in,S2);

% coordinate matrix
M = get_submatrix(coeff1,1:nnz(obj.sparsity_out),I1);

if isa(M,'casadi.Sparsity')
    % coordinates is a dense vector
    res = sum2(M);
    % sparsity of result
    S_res = res;

else
    % apply operator matrix
    res = M*coeff_getnz(S2,coeff2,I2);
    % sparsity of result
    S_res = sparsity(res);
end

% set nonzero coefficients to result
coeffs = coeff_setnz(obj.sparsity_out,find(S_res),res);

% return polynomial pattern
[S,coeffs] = coeff_update(obj.sparsity_out,coeffs);

end

%%
function v = coeff_getnz(S,coeffs,I)
% Return nonzero elements.

nz = sparsity_cast(coeffs,casadi.Sparsity.dense(nnz(S),1));

% select nonzeros
v = nz(I);

end

function coeffs = coeff_setnz(S,I,v)
% Set nonzero elements.

sz = coeff_size(S);
[ii,jj] = coeff_triplet(S);

% new coefficient sparsity
S_coeffs = casadi.Sparsity.triplet(sz(1),sz(2),ii(I),jj(I));

% cast onto nonzeros
coeffs = sparsity_cast(v,S_coeffs);

end
