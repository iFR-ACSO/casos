function [S,coeffs] = coeff_sum(S,coeffs,dim)
% Compute (finish) sum of polynomial coefficient matrix.

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

    % set dimensions of output
    S.matdim(dim) = 1;
    % convert subindices back to linear indices
    ind = sub2ind(size(S),ii,jj);

    coeffs = casadi.Sparsity.triplet(S.nterm,numel(S),it,ind-1);

else
    % finish operation on coefficients
    coeffs = finishMatrixOp(S,coeffs,dim);

    % remove zero terms
    [coeffs,S.degmat,S.indets] = removeZero(coeffs,S.degmat,S.indets);
end

% store coefficients
S = set_coefficients(S,coeffs);

end
