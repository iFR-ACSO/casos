function X = unflatten_sdp(x_flat, n)
    % Convert a flattened symmetric matrix representation back to a full matrix.
    %
    % INPUTS:
    %   x_flat - Vector of independent elements of a symmetric n x n matrix
    %   n      - Dimension of the original matrix
    %
    % OUTPUT:
    %   X      - Reconstructed n x n symmetric matrix

    % Initialize empty symmetric matrix
    X = zeros(n, n);

    % Index counter for flattened vector
    idx = 1;
    
    for i = 1:n
        for j = i:n
            X(i, j) = x_flat(idx);
            X(j, i) = x_flat(idx); % Use symmetry
            idx = idx + 1;
        end
    end
end
