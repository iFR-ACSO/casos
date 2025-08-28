function [S,coeffs] = coeff_dot(obj,S2,coeff1,coeff2)
% Compute coefficient matrix for dot operation.

if ~is_operator(S2)
    % operator evaluation
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

else
    % operator composition
    assert(numel(obj.sparsity_in) == numel(obj.sparsity_out), 'Dimensions mismatch for composition.')

    % find common nonzero coordinates
    [~,I1,I2] = op_intersect(obj.sparsity_in,S2.sparsity_out);

    % coordinate matrices
    M1 = get_submatrix(coeff1,1:nnz(obj.sparsity_out),I1);
    M2 = get_submatrix(coeff2,I2,1:nnz(S2.sparsity_in));

    if isa(M1,'casadi.Sparsity')
        % matrix product on sparsity patterns
        %
        % compute the matrix product
        %
        %   | a11 ... a1L |   | b11 ... b1N |
        %   |  :       :  | * |  :       :  |
        %   | aM1 ... aML |   | bL1 ... bLN |
        %
        % via the element-wise product
        %
        %   | a11 ... a1L |    | b11 ... bL1 |
        %   |  :       :  |    |  :       :  |
        %   | aM1 ... aML |    | b11 ... bL1 |
        %   |  :       :  | .* |  :       :  |
        %   | a11 ... a1L |    | b1N ... bLN |
        %   |  :       :  |    |  :       :  |
        %   | aM1 ... aML |    | b1N ... bLN |
        %   
        A = repmat(M1,size(M2,2),1);
        B = T(reshape(repmat(M2,size(M1,1),1),size(A,2),size(A,1)));

        coeffs = reshape(sum2(A*B),size(M1,1),size(M2,2));

    else
        % composition of operator matrices
        coeffs = M1*M2;
    end

    % remove zero terms
    [coeffs,Si,So] = removeZero(coeffs,S2.sparsity_in,obj.sparsity_out);

    % new operator pattern
    S = new_from_coefficients(coeffs,Si,So);
end

end

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
