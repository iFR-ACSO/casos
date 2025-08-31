function M = get_submatrix(coeffs,I,J)
% Return the IxJ submatrix of the coefficients matrix.

if isa(coeffs,'casadi.Sparsity')
    % subreference of sparsity pattern
    M = sub(coeffs,I-1,J-1);    % CasADi uses 0-index

else
    % subreference coefficient matrix
    M = coeffs(I,J);
end

end
