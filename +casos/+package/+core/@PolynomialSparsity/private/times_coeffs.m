function C = times_coeffs(A,B)
% Multiply coefficient matrices.

if isa(A,'casadi.Sparsity')
    % intersect sparsity
    C = intersect(A,B);

else
    % multiply coefficients element-wise
    C = times(A,B);
end

end
