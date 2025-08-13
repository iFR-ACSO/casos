function [S,coeffs] = coeff_adjoint(obj,coeffs)
% Return coefficient matrix of adjoint.

if is_dual(obj)
    % adjoint of dual
    S = primalize(obj);
    
    % return coefficient matrix for primalization
    coeffs = sparsity_cast(coeffs,coeff_sparsity(S));

else
    % adjoin input-output patterns
    Si = adjoint(obj.sparsity_out);
    So = adjoint(obj.sparsity_in);

    % transpose linear map
    coeffs = T(coeffs);

    % new operator sparsity
    S = new_from_coefficients(coeffs,Si,So);
end
