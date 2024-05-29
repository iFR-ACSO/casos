function c = plus(a,b)
% Add two linear operators.

assert(isequal(size(a),size(b)),'Dimensions of operators mismatch.')

% join input nonzero coordinates
[Sin,I1in,I2in] = op_join(a.sparsity_in,b.sparsity_in);
% join output nonzero coordinates
[Sout,I1out,I2out] = op_join(a.sparsity_out,b.sparsity_out);

% size of resulting operator matrix
sz = [nnz(Sout) nnz(Sin)];

% expand matrices to joint nonzeros
A = expand(a.matrix,sz,I1in,I1out);
B = expand(b.matrix,sz,I2in,I2out);

% new operator
c = a.new_operator(A+B,Sin,Sout);

end
