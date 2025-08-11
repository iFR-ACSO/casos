function coeffs = sparsity_sum(coeffs,sz,dim)
% Compute sum of coefficient sparsity pattern along given dimension.

nt = size(coeffs,1);

% get indices for terms and coefficients
[it,ind] = get_triplet(coeffs);

% convert linear indices to subindices
[ii,jj] = ind2sub(sz,ind+1);

% project coefficients along dimension
switch (dim)
    case 1
        ii = ones(size(jj));
        sz(1) = 1;

    case 2
        jj = ones(size(ii));
        sz(2) = 1;

    otherwise
        error('Invalid dimension input.')
end

% convert subindices back to linear indices
ind = sub2ind(sz,ii,jj);

% build sparsity pattern of coefficient matrix
coeffs = casadi.Sparsity.triplet(nt,prod(sz),it,ind-1);

end
