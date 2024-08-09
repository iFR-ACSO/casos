function [M,Si,So] = op2basis(op,Si,So)
% Return a matrix of nonzero coordinates for given input and output basis
% (sparsity). If no basis is given, the operator's sparsity patterns are
% used.

if nargin < 2
    % return nonzero coordinates
    Si = op.sparsity_in;
    So = op.sparsity_out;
    M = op.matrix;
    return
end

% else
assert(nargin == 3,'Undefined syntax.')
assert(isequal(size(Si),size(op.sparsity_in)),'Input dimensions must not change.')
assert(isequal(size(So),size(op.sparsity_out)),'Output dimensions must not change.')

% join input/output nonzero coordinates
[Si2,Iin] = op_join(op.sparsity_in,Si);
[So2,Iout] = op_join(op.sparsity_out,So);

% expand matrix to joint nonzeros
M0 = expand(op.matrix,[nnz(So2) nnz(Si2)],Iin,Iout);

% find common nonzero coordinates
[~,I1,~] = op_intersect(Si2,Si);
[~,~,I2] = op_intersect(So,So2);

% return nonzeros
M = M0(I2,I1);

end
