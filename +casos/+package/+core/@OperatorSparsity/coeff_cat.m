function [S,coeffs] = coeff_cat(obj,S2,coeff1,coeff2,dim)
% Concatenate operator coefficient matrices.

assert(any(dim==[1 2]), 'Operating dimension must be either 1 or 2.')
assert(is_matrix(obj) && is_matrix(S2), 'Not implemented.')

% detect empty operators
if numel(obj) < 1
    % return second coefficient matrix
    S = S2;
    coeffs = coeff2;
    return

elseif isempty(S2)
    % return first coefficient matrix
    S = casos.Sparsity.create(obj);
    coeffs = coeff1;
    return
end

% else
sza = size(obj);
szb = size(S2);

assert(sza(3-dim) == szb(3-dim),'Dimensions of polynomials being concatenated are not consistent.')

switch (dim)
    case 1
        % vertical concatenation
        %
        %   | A |
        %   | B |
        %
        % corresponds to the operator
        %
        %   (x) -> (Ax, Bx)
        %

        % join input nonzero coordinates
        [Si,I1,I2] = op_join(obj.sparsity_in,S2.sparsity_in);
        
        % expand coefficient matrices
        coeffA = expand_matrix(obj.sparsity_M,coeff1,[nnz(obj.sparsity_out) nnz(Si)],I1,1:nnz(obj.sparsity_out));
        coeffB = expand_matrix(S2.sparsity_M,coeff2,[nnz(S2.sparsity_out) nnz(Si)],I2,1:nnz(S2.sparsity_out));
        
        % concatenate output sparsity patterns
        So = vertcat(obj.sparsity_out,S2.sparsity_out);

        % vertically concatenate coefficient matrices
        coeffs = vertcat(coeffA,coeffB);

    case 2
        % horizontal concatenation
        %
        %   | A  B |
        %
        % corresponds to the operator
        %
        %   (x,y) -> (Ax + By)
        %

        % join output nonzero coordinates
        [So,I1,I2] = op_join(obj.sparsity_out,S2.sparsity_out);

        % expand coefficient matrices
        coeffA = expand_matrix(obj.sparsity_M,coeff1,[nnz(So) nnz(obj.sparsity_in)],1:nnz(obj.sparsity_in),I1);
        coeffB = expand_matrix(S2.sparsity_M,coeff2,[nnz(So) nnz(S2.sparsity_in)],1:nnz(S2.sparsity_in),I2);

        % concatenate input sparsity patterns
        Si = vertcat(obj.sparsity_in,S2.sparsity_in);

        % horizontally concatenate coefficient matrices
        coeffs = horzcat(coeffA,coeffB);
end

% remove zero terms
[coeffs,Si,So] = removeZero(coeffs,Si,So);

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
