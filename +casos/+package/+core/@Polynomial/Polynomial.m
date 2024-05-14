classdef (Abstract) Polynomial < casos.package.core.GenericPolynomial
% Polynomials with constant and symbolic coefficients.

properties (Access=protected)
    % coefficient matrix
    % see casos.Sparsity for details
    coeffs = [];
end

methods (Static, Abstract, Access=protected)
    c = new_coeff(varargin);
    s = sym_coeff(varargin);

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
            assert(nargin > 1,'Not enough arguments.')
            assert(nargin < 3,'Too many arguments.')

            psparsity = casos.Sparsity(varargin{1});
            coeffs = obj.new_coeff(coeff_sparsity(psparsity),varargin{2});

        elseif ischar(varargin{1}) || isa(varargin{1},'casos.Indeterminates')
            % indeterminate (pvar / mpvar syntax)
            [vars,I] = sort(casos.Indeterminates(varargin{:}));
            psparsity = to_vector(casos.Sparsity(vars),I);
            coeffs = obj.new_coeff(coeff_sparsity(psparsity),1);

        else
            % zero-degree polynomial (casadi syntax)
            C = obj.new_coeff(varargin{:});
            [psparsity,coeffs] = casos.Sparsity.coeff_zerodegree(C);
        end

        obj = set_sparsity(obj,psparsity);
        obj.coeffs = coeffs;
    end

    %% Getter
    function tf = is_symbolic(obj)
        % Check if polynomial has symbolic coefficients.
        tf = is_symbolic(obj.coeffs);
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

    function [tf,I] = is_monom(obj)
        % Check if polynomial is vector of monomials (legacy).
        [tf0,I] = is_monom(obj.get_sparsity);
        tf = (tf0 && is_coeff_one(obj));
    end

    function [tf,I] = is_indet(obj)
        % Check if polynomial is a vector of indeterminates.
        [tf0,I] = is_monom(obj);
        tf = (tf0 && is_homogeneous(obj,1));
    end
end

methods
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

    function v = casos.Indeterminates(obj)
        % Convert to indeterminates.
        [tf,I] = is_indet(obj);

        assert(tf,'Polynomial must be a vector of indeterminates.')

        % get indeterminate variables
        vars = indeterminates(obj);
        % return variables in order
        v = vars(I);
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
        % Display polynomial as matrix.
        disp_matrix(obj);
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
