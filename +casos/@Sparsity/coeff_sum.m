function [S,coeffs] = coeff_sum(S,coeffs,dim)
% Compute sum of polynomial coefficient matrix.

sz = sizeofMatrixOp(S,dim);

if isa(coeffs,'casadi.Sparsity')
    % compute sparsity pattern
    coeffs = sparsity_sum(coeffs,size(S),dim);

else
    % prepare for operation on coefficients
    cfa = prepareMatrixOp(S,coeffs,dim);

    % sum coefficients
    cfb = evalMatrixOp(cfa,@sum,dim);
    
    % finish operation on coefficients
    coeffs = finishMatrixOp(cfb,sz,dim);

    % remove zero terms
    [coeffs,S.degmat,S.indets] = removeZero(coeffs,S.degmat,S.indets);
end

% set dimensions of output
S.matdim = sz;
% store coefficients
S = set_coefficients(S,coeffs);

end
