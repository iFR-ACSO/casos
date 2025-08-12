classdef (Abstract) Polynomial < casos.package.core.GenericPolynomial
% Polynomials with constant and symbolic coefficients.

properties (Access=protected)
    % coefficient matrix
    % see casos.Sparsity for details
    coeffs = [];

    % polynomial sparsity pattern
    poly_sparsity;
end

methods (Static, Abstract, Access=protected)
    c = new_coeff(varargin);
    p = new_poly(varargin);
end

methods (Static, Abstract)
    % Static constructors
    p = sym(dstr,varargin);
    p = empty();
    p = zeros(varargin);
    p = ones(varargin);
    p = eye(varargin);
end

methods (Abstract, Access=protected)
    tf = is_coeff_one(obj);
end

methods
    %% Public constructor
    function obj = Polynomial(varargin)
        % Create polynomial variable.
        if nargin == 0
            % empty polynomial
            psparsity = casos.Sparsity(0,0);
            coeffs = obj.new_coeff;

        elseif isa(varargin{1},'casos.package.core.Polynomial')
            % copy or convert polynomial
            assert(nargin == 1,'Too many arguments.')

            psparsity = sparsity(varargin{1});
            coeffs = obj.new_coeff(varargin{1}.coeffs);

        elseif isa(varargin{1},'casos.Sparsity')
            % sparsity pattern syntax
            assert(nargin < 3,'Too many arguments.')

            psparsity = casos.Sparsity(varargin{1});
            coeffs = obj.new_coeff(coeff_sparsity(psparsity),varargin{2:end});

        elseif ischar(varargin{1}) || isa(varargin{1},'casos.Indeterminates')
            % indeterminate (pvar / mpvar syntax)
            [vars,~,I] = sort(casos.Indeterminates(varargin{:}));
            psparsity = to_vector(casos.Sparsity(vars),I,isrow(vars));
            coeffs = obj.new_coeff(coeff_sparsity(psparsity),1);

        else
            % zero-degree polynomial (casadi syntax)
            % or operator syntax
            C = obj.new_coeff(varargin{1});
            psparsity = casos.Sparsity(sparsity(C),varargin{2:end});
            coeffs = sparsity_cast(C, coeff_sparsity(psparsity));
        end

        obj = set_sparsity(obj,psparsity);
        obj.coeffs = coeffs;
    end

    %% Getter
    function tf = is_symbolic(obj)
        % Check if polynomial has symbolic coefficients.
        tf = is_symbolic(obj.coeffs(coeff_find(obj.get_sparsity)));
    end

    function tf = is_symexpr(obj)
        % Check if polynomial contains symbolic expressions.
        tf = ~is_constant(obj.coeffs);
    end

    function tf = is_zero(obj)
        % Check if polynomial is equal to zero.
        tf = (is_zerodegree(obj) && is_zero(obj.coeffs));
    end

    function tf = is_one(obj)
        % Check if polynomial is equal to identity.
        tf = (is_zerodegree(obj) && is_coeff_one(obj));
    end

    function tf = is_constant(obj)
        % Check if polynomial is constant.
        tf = (is_zerodegree(obj) && ~is_symexpr(obj));
    end

    function tf = is_monom(obj)
        % Check if polynomial is vector of monomials (legacy).
        tf = (is_monom(obj.get_sparsity) && is_coeff_one(obj));
    end

    function tf = is_indet(obj)
        % Check if polynomial is a vector of indeterminates.
        tf = (is_monom(obj) && is_homogeneous(obj,1));
    end

    function tf = is_equal(obj,p)
        % Check if polynomials are equal.
        tf = is_equal@casos.package.core.GenericPolynomial(obj,p) ...
            && is_equal(obj.coeffs,p.coeffs);
    end

    function l = list_of_coeffs(obj)
        % Return a list of coefficients.
        l = coeff_list(obj.get_sparsity,obj.coeffs);
    end
end

methods
    % check well-posedness
    tf = is_wellposed(obj);

    % public RedefinesParen interface
    p = cat(dim,varargin);

    function obj = reshape(obj,varargin)
        % Reshape polynomial matrix.
        assert(length(varargin{1}) <= 2, 'Size vector must not exceed two elements.')
        assert(length(varargin) <= 2, 'Size arguments must not exceed two scalars.')

        S = reshape(obj.get_sparsity,varargin{:});
        obj = obj.set_sparsity(S);
    end

    function obj = repmat(obj,varargin)
        % Repeat polynomial matrix.
        [S,obj.coeffs] = coeff_repmat(obj.get_sparsity,obj.coeffs,varargin{:});
        obj = obj.set_sparsity(S);
    end

    function obj = project(obj,S)
        % Project onto sparsity pattern.
        if ~check_sz_equal(obj,S)
            throw(casos.package.core.IncompatibleSizesError.other(obj,S));
        end

        [S,obj.coeffs] = coeff_project(obj.get_sparsity,obj.coeffs,casos.Sparsity(S));
        obj = obj.set_sparsity(S);
    end

    function obj = sparsity_cast(obj,S)
        % Cast sparsity pattern onto polynomial.
        assert(nnz(obj) == nnz(S), 'Mismatching number of nonzero coefficients.')

        S = casos.Sparsity(S);

        obj.coeffs = sparsity_cast(obj.coeffs,coeff_sparsity(S));
        obj = obj.set_sparsity(S);
    end

    function v = casos.Indeterminates(obj)
        % Convert to indeterminates.
        assert(is_indet(obj),'Polynomial must be a vector of indeterminates.')

        v = vector_to_indeterminates(obj.poly_sparsity);
    end

    function M = full(obj)
        % Convert to full matrix.
        assert(is_zerodegree(obj), 'Can only convert polynomial of degree zero.')

        M = full(reshape(obj.coeffs,size(obj)));
    end

    %% Unary operators
    function p = uplus(p)
        % Unary plus.
        p.coeffs = uplus(p.coeffs);
    end

    function p = uminus(p)
        % Unary minus.
        p.coeffs = uminus(p.coeffs);
    end

    function c = minus(a,b)
        % Substract two polynomials.
        c = plus(a, uminus(b));
    end

    %% Display
    function disp(obj)
        % Display polynomial.
        if ~is_operator(obj)
            % display as matrix
            disp_matrix(obj);
        else
            % display operator
            disp@casos.package.core.GenericPolynomial(obj);
        end
    end
end

methods (Access=protected)
    % helper for static constructors
    obj = new_sym(obj,dstr,varargin);

    % protected RedefinesParen interface
    obj = parenAssign(obj,idx,varargin);
    % obj = parenDelete(obj,idx);
    % n = parenListLength(obj,idx,context);
    varargout = parenReference(obj,index);
end

end
