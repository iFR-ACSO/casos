function [S,coeffs] = coeff_sum(S,coeffs,dim)
% Compute sum of polynomial coefficient matrix.

sz = sizeofMatrixOp(S,dim);

if isa(coeffs,'casadi.Sparsity')
    % compute sparsity pattern
    [it,ind] = get_triplet(coeffs);
    % convert linear indices to subindices
    [ii,jj] = ind2sub(size(S),ind+1);

    % project coefficients along dimension
    switch (dim)
        case 1
            ii = ones(size(jj));

        case 2
            jj = ones(size(ii));

        otherwise
            error('Invalid dimension input.')
    end

    % convert subindices back to linear indices
    ind = sub2ind(sz,ii,jj);

    coeffs = casadi.Sparsity.triplet(S.nterm,prod(sz),it,ind-1);

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
