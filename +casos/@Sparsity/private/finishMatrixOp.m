function coeffs = finishMatrixOp(coeffs,sz,dim)
% Finish after matrix operation along specified dimension.
%
% This function reshapes the coefficient matrix C after applying
% matrixop(C,dim) to store in PS object.

% The coefficients of a MxN matrix polynomial are stored as
%
%       | a11 ... aM1 ... a1N ... aMN |
%   D = |  :       :       :       :  |
%       | z11 ... zM1 ... z1N ... zMN |
%
% where each row of D corresponds to a term of the polynomial.

if isequal(dim,'all')
    % Matrix operation over all elements (scalar)
    % nothing to do
    return
end

% else:
ne = size(coeffs,3-dim);

switch (dim)
    case 1
        % Matrix operation along first dimension (row vector)
        % input is
        %
        %   D = | a1 ... aN ... z1 ... zN |
        coeffs = T(reshape(coeffs,prod(sz),ne/sz(2)));

    case 2
        % Matrix operation along second dimension (column vector)
        % input is
        %
        %       | a1 |
        %       |  : |
        %       | z1 |
        %   D = |  : |
        %       | aM |
        %       |  : |
        %       | zM |
        %
        coeffs = reshape(coeffs,ne/sz(1),prod(sz));

    otherwise
        error('Invalid dimension input.')
end

end
