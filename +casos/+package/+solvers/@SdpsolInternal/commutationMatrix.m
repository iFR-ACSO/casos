function K = commutationMatrix(~, n)
% commutationMatrix: Construct commutation matrix K_{n,n}.
%   K = commutationMatrix(n) returns the sparse commutation matrix of size
%   n^2-by-n^2 such that:
%       K * vec(X) = vec(X.')

% Preallocate indices
rows = zeros(n^2,1);
cols = zeros(n^2,1);

% Compute permutation mapping
idx = 1;
for i = 1:n
    for j = 1:n
        % vec(X) index for (i,j)
        col = (j-1)*n + i;
        % vec(X.') index for (i,j)
        row = (i-1)*n + j;

        rows(idx) = row;
        cols(idx) = col;
        idx = idx + 1;
    end
end

% Build sparse permutation matrix
K = sparse(rows, cols, 1, n^2, n^2);
end