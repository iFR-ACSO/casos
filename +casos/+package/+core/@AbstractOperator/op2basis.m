function [M,S] = op2basis(op,Si,So)
% Return a matrix of nonzero coordinates for given input and output basis
% (sparsity). If no basis is given, the operator's sparsity patterns are
% used.

if nargin < 2
    % return nonzero coordinates
    S = sparsity(op);
    M = op.matrix;
    return

elseif nargin < 3
    % OperatorSparsity syntax
    S = Si;
    assert(isa(S,'casos.package.core.OperatorSparsity'),'Second input must be operator sparsity pattern.')

    % get operator for input-output patterns
    M0 = op2basis(op,S.sparsity_in,S.sparsity_out);
    % project onto linear map pattern
    M = project(M0,matrix_sparsity(S));
    return
end

% else
assert(isequal(size(Si),size_in(op)),'Input dimensions must not change.')
assert(isequal(size(So),size_out(op)),'Output dimensions must not change.')

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

% return operator sparsity pattern
S = casos.package.core.OperatorSparsity(sparsity(M),Si,So);

end
