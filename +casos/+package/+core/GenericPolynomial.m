classdef (Abstract) GenericPolynomial < casos.package.core.AlgebraicObject ...
        & casos.package.core.PolynomialInterface
% Base class for all polynomial objects.

properties (Access=private)
    poly_sparsity;   % polynomial sparsity pattern
end

properties (Dependent)
    nvars;
    nterm;
    maxdeg;
    mindeg;
end

methods
    %% Superclass constructor
    function obj = GenericPolynomial(sparsity)
        % Create a new object with sparsity pattern.
        if nargin > 0
            obj.poly_sparsity = sparsity;
        else
            obj.poly_sparsity = casos.Sparsity;
        end
    end

    %% Getter
    function n = get.nvars(obj)
        % Number of indeterminate variables.
        n = obj.poly_sparsity.nvars;
    end

    function n = get.nterm(obj)
        % Number of monomials.
        n = obj.poly_sparsity.nterm;
    end

    function d = get.mindeg(obj)
        % Minimum degree of polynomial.
        d = obj.poly_sparsity.mindeg;
    end

    function d = get.maxdeg(obj)
        % Maximum degree of polynomial.
        d = obj.poly_sparsity.maxdeg;
    end

    function n = nnz(obj)
        % Return number of nonzero coefficients.
        n = nnz(obj.poly_sparsity);
    end

    function n = numel(obj)
        % Return number of elements.
        n = numel(obj.poly_sparsity);
    end

    function varargout = size(obj,varargin)
        % Return size of polynomial.
        [varargout{1:nargout}] = size(obj.poly_sparsity,varargin{:});
    end

    function x = indeterminates(obj)
        % Return indeterminate variables of polynomial.
        x = casos.Indeterminates(obj.poly_sparsity);
    end

    function tf = isrow(obj)
        % Check if polynomial is a row vector.
        tf = isrow(obj.poly_sparsity);
    end

    function tf = iscolumn(obj)
        % Check if polynomial is a column vector.
        tf = iscolumn(obj.poly_sparsity);
    end

    function tf = isvector(obj)
        % Check if polynomial is a vector.
        tf = isvector(obj.poly_sparsity);
    end

    function tf = isscalar(obj)
        % Check if polyonomial is a scalar.
        tf = isscalar(obj.poly_sparsity);
    end

    function tf = is_zerodegree(obj)
        % Check if polynomial is of degree zero.
        tf = is_zerodegree(obj.poly_sparsity);
    end

    function S = sparsity(obj)
        % Return (copy of) sparsity pattern.
        S = casos.Sparsity(obj.poly_sparsity);
    end
end

methods (Access=protected)
    function obj = set_sparsity(obj,sparsity)
        % Set polynomial sparsity pattern.
        obj.poly_sparsity = sparsity;
    end

    function sparsity = get_sparsity(obj)
        % Return polynomial sparsity pattern.
        sparsity = obj.poly_sparsity;
    end
end

end
