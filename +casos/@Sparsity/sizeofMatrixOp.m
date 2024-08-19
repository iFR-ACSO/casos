function sz = sizeofMatrixOp(S,dim)
% Return size of matrix operation along specified dimension.

sz = S.size;

switch (dim)
    case 'all'
        % Matrix operation over all elements
        % result is scalar
        sz = [1 1];

    case 1
        % Matrix operation along first dimension
        % result is row vector
        sz(1) = 1;

    case 2
        % Matrix operation along second dimension
        % result is column vector
        sz(2) = 1;

    otherwise
        error('Invalid dimension input.')
end

end
