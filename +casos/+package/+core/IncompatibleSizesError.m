% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef IncompatibleSizesError < MException
% Error thrown if inputs have compatible sizes for an operation.

methods (Access=protected)
    function err = IncompatibleSizesError(id,msg,a,b)
        % New error with identifier.
        err@MException(id,msg,size2str(a),size2str(b));
    end
end

methods
    function msgtext = getReport(err,varargin)
        % Return formatted error message text.
        if nargin < 2, varargin = {'basic'}; end % omit stack from report

        msgtext = getReport@MException(err,varargin{:});
    end
end

methods (Static)
    function err = assign(a,b)
        % New error for assignment operation.
        err = casos.package.core.IncompatibleSizesError('MATLAB:subsassigndimmismatch', ...
                ['Unable to perform assignment because the size of the left polynomial is [%s]' ...
                ' and the size of the right polynomial is [%s].'],a,b);
    end

    function err = basic(a,b)
        % New error for basic operation.
        err = casos.package.core.IncompatibleSizesError('MATLAB:sizeDimensionsMustMatch', ...
                'Polynomials have incompatible sizes for this operation: [%s] vs. [%s].',a,b);
    end

    function err = concat(a,b)
        % New error for concatenation.
        err = casos.package.core.IncompatibleSizesError('MATLAB:catenate:dimensionMismatch', ...
                'Dimensions of polynomials being concatenated are not consistent: [%s] vs. [%s].',a,b);
    end

    function err = matrix(a,b)
        % New error for matrix operation.
        err = casos.package.core.IncompatibleSizesError('MATLAB:innerdim', ...
                'Incorrect dimensions for polynomial matrix multiplication: [%s] vs. [%s].',a,b);
    end

    function err = mpower(a,b)
        % New error for matrix power.
        err = casos.package.core.IncompatibleSizesError('MATLAB:mpower:notScalarAndSquareMatrix', ...
                'Incorrect dimensions for raising a polynomial matrix to a power: [%s] and [%s].',a,b);
    end

    function err = substitute(a,b)
        % New error for substitution.
        err = casos.package.core.IncompatibleSizesError('', ...
                ['Unable to perform substitution because the size of the indeterminate variables is [%s]' ...
                ' and the size of the polynomial to be substituted into is [%s].'],a,b);
    end

    function err = other(a,b)
        % New error for other operations.
        err = casos.package.core.IncompatibleSizesError('', ...
                'Polynomial dimension mismatch: [%s] vs. [%s].',a,b);
    end
end

end

function str = size2str(obj)
    % Convert array size to string.
    str = sprintf('%dx%d', size(obj,1), size(obj,2));
end
