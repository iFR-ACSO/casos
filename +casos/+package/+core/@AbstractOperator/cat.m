function c = cat(dim,a,b)
% Concatenate linear operators.

assert(iscolumn(a) && iscolumn(b), 'Not implemented.')

switch (dim)
    case 0
        % Diagonal concatenation
        assert(isequal(size_in(a,2),size_in(b,2)),'Input dimensions mismatch.')
        assert(isequal(size_out(a,2),size_out(b,2)),'Output dimensions mismatch.')

        Sin = vertcat(a.sparsity_in,b.sparsity_in);
        Sout = vertcat(a.sparsity_out,b.sparsity_out);

        % block concatenation
        M = blkdiag(a.matrix,b.matrix);

    case 1
        % Vertical concatenation
        % y = [A*x; B*x]
        assert(isequal(size_in(a),size_in(b)),'Input dimensions mismatch.')
        assert(isequal(size_out(a,2),size_out(b,2)),'Output dimensions mismatch.')

        % join input nonzero coordinates
        [Sin,I1in,I2in] = op_join(a.sparsity_in,b.sparsity_in);
        Sout = vertcat(a.sparsity_out,b.sparsity_out);

        % expand matrices to joint nonzeros
        A = expand(a.matrix,[nnz_out(a) nnz(Sin)],I1in,1:nnz_out(a));
        B = expand(b.matrix,[nnz_out(b) nnz(Sin)],I2in,1:nnz_out(b));

        % y = [A; B]*x
        M = vertcat(A,B);

    case 2
        % Horizontal concatenation
        % y = A*x1 + B*x2
        assert(isequal(size_in(a,2),size_in(b,2)),'Input dimensions mismatch.')
        assert(isequal(size_out(a),size_out(b)),'Output dimensions mismatch.')

        % join output nonzero coordinates
        [Sout,I1out,I2out] = op_join(a.sparsity_out,b.sparsity_out);
        Sin = vertcat(a.sparsity_in,b.sparsity_in);

        % expand matrices to joint nonzeros
        A = expand(a.matrix,[nnz(Sout) nnz_in(a)],1:nnz_in(a),I1out);
        B = expand(b.matrix,[nnz(Sout) nnz_in(b)],1:nnz_in(b),I2out);

        % y = [A B]*[x1;x2]
        M = horzcat(A,B);

end
   
% new operator
c = a.new_operator(M,Sin,Sout);

end
