function [coeffs,sz] = prepareMatrixOp(coeffs,sz,dim)
% Prepare for matrix operation along specified dimension.
%
% This function returns the coefficient matrix C of a matrix polynomial in
% such a form that matrixop(C,dim) corresponds to the same matrix operation
% for each coefficient of the polynomial.

% The coefficients of a MxN matrix polynomial are stored as
%
%       | a11 ... aM1 ... a1N ... aMN |
%   C = |  :       :       :       :  |
%       | z11 ... zM1 ... z1N ... zMN |
%
% where each row of C corresponds to a term of the polynomial.
nt = size(coeffs,1);

switch (dim)
    case 'all'
        % Matrix operation over all elements; return C 
        % result is scalar
        sz = [1 1];

    case 1
        % Matrix operation along first dimension; return
        %
        %       | a11 ... a1N ... z11 ... z1N |
        %   D = |  :       :       :       :  |
        %       | aM1 ... aMN ... zM1 ... zMN |
        %
        coeffs = reshape(coeffs',sz(1),sz(2)*nt);
        
        % result is row vector
        sz(1) = 1;

    case 2
        % Matrix operation along second dimension; return
        %
        %       | a11 ... a1N |
        %       |  :       :  |
        %       | z11 ... z1N |
        %   D = |  :       :  |
        %       | aM1 ... aMN |
        %       |  :       :  |
        %       | zM1 ... zMN |
        %
        coeffs = reshape(coeffs,sz(1)*nt,sz(2));

        % result is column vector
        sz(2) = 1;

    otherwise
        error('Invalid dimension input.')
end

end
