classdef (Abstract) PolynomialInterface < casos.package.core.Printable
% Base class for polynomial-like objects.

methods (Abstract)
    tf = is_wellposed(obj);
end

methods
    %% Concatenation interface
    function p = horzcat(varargin)
        % Horizontal concatenation.
        p = cat(2,varargin{:});
    end

    function p = vertcat(varargin)
        % Vertical concatenation.
        p = cat(1,varargin{:});
    end

    function p = blkdiag(varargin)
        % Block diagonal concatenation.
        p = cat(0,varargin{:});
    end

    function p = cat(dim,varargin)
        % Generic concatenation.
        switch (nargin-1)
            case 0, p = [];
            case 1, p = varargin{:};
            case 2, error('Notify the developers.')
            otherwise
                % more than two inputs
                N = (nargin-1)/2;
                % recursion
                p1 = cat(dim,varargin{1:floor(N)});
                p2 = cat(dim,varargin{floor(N)+1:end}); 
                % concatenate
                p = cat(dim,p1,p2);
        end
    end
end

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

    function tf = check_sz_assign(a,b)
        % Check if sizes are compatible for matrix assignment.
        sza = size(a);
        szb = size(b);

        % dimensions are compatible if equal or right side is row/column
        tf = all(sza == szb | szb == 1);
    end

    function tf = check_sz_mtimes(a,b)
        % Check if sizes are compatible for matrix multiplication.
        tf = size(a,2) == size(b,1);
    end
end

end
