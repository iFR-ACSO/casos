classdef (Abstract) PolynomialInterface < casos.package.core.Printable
% Base class for polynomial-like objects.

methods (Access=protected)
    %% Utilities
    function tf = check_sz_equal(a,b)
        % Check if two objects have same size.
        tf = isequal(size(a),size(b));
    end

    function tf = check_sz_comptbl(a,b)
        % Check if sizes are compatible for basic (array) operations.
        % input dimensions
        sza = size(a);
        szb = size(b);
        
        % compare dimensions
        I = (sza == szb);
        % find one dimension
        I1 = (sza == 1) | (szb == 1);
        
        % dimensions are compatible if equal or one summand is row/column
        tf = all(I | I1);
    end

    function tf = check_sz_mtimes(a,b)
        % Check if sizes are compatible for matrix multiplication.
        tf = size(a,2) == size(b,1);
    end
end

end
