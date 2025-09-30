classdef IncompatibleSizesError < MException
% Error thrown if inputs have compatible sizes for an operation.

methods (Access=protected)
    function err = IncompatibleSizesError(id,a,b)
        % New error with identifier.
        errsz = 'Polynomials have incompatible sizes for this operation ([%s] vs. [%s]).';
        
        err@MException(id,errsz,size2str(a),size2str(b));
    end
end

methods (Static)
    function err = basic(a,b)
        % New error for basic operation.
        err = casos.package.core.IncompatibleSizesError('MATLAB:sizeDimensionsMustMatch',a,b);
    end

    function err = matrix(a,b)
        % New error for matrix operation.
        err = casos.package.core.IncompatibleSizesError('MATLAB:innerdim',a,b);
    end

    function err = mdivide(a,b)
        % New error for matrix division.
        err = casos.package.core.IncompatibleSizesError('MATLAB:dimagree',a,b);
    end

    function err = other(a,b)
        % New error for other operations.
        err = casos.package.core.IncompatibleSizesError('',a,b);
    end
end

end

function str = size2str(obj)
    % Convert array size to string.
    str = sprintf('%dx%d', size(obj,1), size(obj,2));
end
