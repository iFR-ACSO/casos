classdef (Abstract) GenericPolynomial ...
        < casos.package.core.AlgebraicObject & casos.package.core.PolynomialInterface
% Base class for all polynomial objects.

properties (Abstract, Access=protected)
    poly_sparsity;   % polynomial sparsity pattern
end

properties (Dependent)
    nvars;
    nterm;
    maxdeg;
    mindeg;
end

methods
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

    function varargout = degrees(obj)
        % Return degrees of polynomial.
        [varargout{1:nargout}] = degrees(obj.poly_sparsity);
    end

    function x = indeterminates(obj)
        % Return indeterminate variables of polynomial.
        x = indeterminates(obj.poly_sparsity);
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

    function tf = is_homogeneous(obj,varargin)
        % Check if polynomial is homogeneous.
        tf = is_homogeneous(obj.poly_sparsity,varargin{:});
    end

    function tf = is_dense(obj)
        % Check if polynomial has no sparse coefficients.
        tf = is_dense(obj.poly_sparsity);
    end

    function tf = is_equal(obj,p)
        % Check if polynomials are equal.
        tf = is_equal(obj.poly_sparsity,p.poly_sparsity);
    end

    function S = sparsity(obj)
        % Return (copy of) sparsity pattern.
        S = casos.Sparsity(obj.poly_sparsity);
    end

    function Z = grambasis(obj)
        % Return a Gram basis for this polynomial.
        Z = grambasis(obj.poly_sparsity);
    end

    function l = list_of_degree(obj)
        % Return a list of degrees.
        l = list_of_degree(obj.poly_sparsity);
    end

    function l = list_of_indets(obj)
        % Return a list of indeterminate variables.
        l = list_of_indets(obj.poly_sparsity);
    end

    function Z = monomials(obj,deg)
        % Return monomial sparsity pattern.
        if nargin > 1
            % create pattern
            Z = monomials@casos.package.core.AlgebraicObject(obj,deg);
        else
            % return pattern
            Z = monomials(obj.poly_sparsity);
        end
    end

    function tf = is_wellposed(obj)
        % Check if polynomial is well posed.
        tf = is_wellposed(obj.poly_sparsity);
    end

    %% Conversion
    function v = casos.Indeterminates(obj) %#ok<STOUT,MANU>
        % Convert to indeterminates.
        error('Notify the developers.')
    end

    function f = to_sxfunction(obj,varargin)
        % Return casadi.Function object using SX.
        X = casadi.SX.sym('x',obj.nvars,1);
        p = subs(obj,indeterminates(obj),casos.package.polynomial(X));
        f = casadi.Function('f',num2cell(X),{casadi.SX(p)},str(obj.indeterminates),{'poly'},varargin{:});
    end

    function f = to_mxfunction(obj,varargin)
        % Return casadi.Function object using MX.
        X = casadi.MX.sym('x',obj.nvars,1);
        p = subs(obj,indeterminates(obj),casos.package.polynomial(X));
        f = casadi.Function('f',num2cell(X),{casadi.MX(p)},str(obj.indeterminates),{'poly'},varargin{:});
    end

    function f = to_function(obj,varargin)
        % Return casadi.Function object.
        % Overload for default behaviour.
        f = to_sxfunction(obj,varargin{:});
    end

    %% Concatenation
    function p = cat(dim,varargin)
        % Generic concatenation.
        p = cat@casos.package.core.PolynomialInterface(dim,varargin{:});
    end

    function p = vertcat(varargin)
        % Vertical concatenation.
        p = vertcat@casos.package.core.PolynomialInterface(varargin{:});
    end

    function p = horzcat(varargin)
        % Horizontal concatenation.
        p = horzcat@casos.package.core.PolynomialInterface(varargin{:});
    end
end

methods (Access=protected)
    %% Protected interface
    function obj = set_sparsity(obj,sparsity)
        % Set polynomial sparsity pattern.
        obj.poly_sparsity = sparsity;
    end

    function sparsity = get_sparsity(obj)
        % Return polynomial sparsity pattern.
        sparsity = obj.poly_sparsity;
    end

    function varargout = parenReference(obj,indexOp)
        % Handle getters on referenced polynomial.
        I = logical(sparse(size(obj,1),size(obj,2)));
        I.(indexOp(1)) = true;

        assert(length(indexOp) > 1, 'Notify the developers.')

        if indexOp(2).Type == "Dot"
            name = indexOp(2).Name;
        else
            name = '';
        end

        % handle dot reference
        switch (name)
            case 'mindeg'
                res = min(get_degree(obj.poly_sparsity,I));
            case 'maxdeg'
                res = max(get_degree(obj.poly_sparsity,I));
            case 'nvars'
                res = get_nvars(obj.poly_sparsity,I);
            case 'nterm'
                res = get_nterm(obj.poly_sparsity,I);
            case 'degrees'
                if length(indexOp) > 2, res = degrees(obj.poly_sparsity,I);
                else, [varargout{1:nargout}] = degrees(obj.poly_sparsity,I);
                    return
                end
            case 'indeterminates'
                res = get_indets(obj.poly_sparsity,I);
            case 'monomials'
                res = get_monoms(obj.poly_sparsity,I);
            otherwise
                % getter not supported
                p = obj.(indexOp(1));
                [varargout{1:nargout}] = p.(indexOp(2:end));
                return
        end

        if length(indexOp) > 2
            [varargout{1:nargout}] = res.(indexOp(3:end));
        else
            varargout = {res};
        end
    end

    function obj = parenAssign(obj,indexOp,varargin)
        % Handle forwarded assignments to subreference.
        idx = indexOp(1);
        
        if length(indexOp) > 1
            % perform parantheses reference
            p = obj.(idx);
            % forward assign
            [p.(indexOp(2:end))] = varargin{:};
            % re-assign modified element
            obj.(idx) = p;
        else
            error('Notify the developers.')
        end
    end
end

methods (Static,Access=protected)
    %% Static protected interface
    function S = sym_pattern(w,sz,type)
        % Return a sparsity pattern for symbolic polynomials.
        if nargin < 1
            % default syntax
            S = casos.Sparsity.scalar;
        elseif isnumeric(w)
            % constant polynomial
            assert(nargin < 2, 'Undefined syntax.')
            S = casos.Sparsity.dense(w);
        elseif nargin < 2
            % sparsity pattern
            S = casos.Sparsity(w);
        elseif ischar(sz)
            % type given
            assert(nargin < 3, 'Undefined syntax.')
            assert(isequal(sz,'gram'), 'Type "%s" undefined.', sz)
            S = gram(casos.Sparsity(w));
        elseif nargin < 3
            % monomials and size given
            S = casos.Sparsity.dense(sz,w);
        else
            assert(nargin < 4, 'Undefined syntax.')
            assert(isequal(type,'gram'), 'Type "%s" undefined.', sz)
            S = gram(casos.Sparsity.dense(sz,w));
        end
    end
end

end
