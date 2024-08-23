function [M,S] = op2basis(op,S)
% Return a matrix of nonzero coordinates for given input and output basis
% (sparsity). If no basis is given, the operator's sparsity patterns are
% used.

if nargin < 2
    % return nonzero coordinates
    S = sparsity(op);
    M = op.matrix;
    return
end

% else
assert(isequal(size_in(S),size_in(op)),'Input dimensions must not change.')
assert(isequal(size_out(S),size_out(op)),'Output dimensions must not change.')

% join input/output nonzero coordinates
[Si2,Iin] = op_join(op.sparsity_in,S.sparsity_in);
[So2,Iout] = op_join(op.sparsity_out,S.sparsity_out);

% expand matrix to joint nonzeros
M0 = expand(op.matrix,[nnz(So2) nnz(Si2)],Iin,Iout);

% find common nonzero coordinates
[~,I1,~] = op_intersect(Si2,S.sparsity_in);
[~,~,I2] = op_intersect(S.sparsity_out,So2);

% return nonzeros
M = project(M0(I2,I1),matrix_sparsity(S));

end
