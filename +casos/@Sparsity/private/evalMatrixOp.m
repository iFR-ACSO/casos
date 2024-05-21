function coeffs = evalMatrixOp(coeffs,op,dim)
% Evaluate matrix operation along specified dimension.

switch (dim)
    case 'all'
        % matrix operation over all elements
        coeffs = op(coeffs,2);

    case {1 2}
        % matrix operation along dimension
        coeffs = op(coeffs,dim);

    otherwise
        error('Invalid dimension input.')
end

end
