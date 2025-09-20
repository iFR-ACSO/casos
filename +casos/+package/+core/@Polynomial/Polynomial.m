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
    p = eye(varargin);
    p = ones(varargin);
    p = zeros(varargin);
    p = id_operator(varargin);
    p = one_operator(varargin);
    p = zero_operator(varargin);
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

        elseif isa(varargin{1},'casos.package.core.Polynomial') && nargin < 2
            % copy or convert polynomial
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

        elseif nargin == 3 && isa(varargin{3},'casos.Sparsity')
            % operator syntax
            C = obj.new_coeff(varargin{1});
            psparsity = casos.Sparsity(sparsity(C),varargin{2:end});
            coeffs = sparsity_cast(C, coeff_sparsity(psparsity));

        else
            % zero-degree polynomial (casadi syntax)
            C = obj.new_coeff(varargin{:});
            psparsity = casos.Sparsity(sparsity(C));
            coeffs = sparsity_cast(C, coeff_sparsity(psparsity));
        end

        obj = set_sparsity(obj,psparsity);
        obj.coeffs = coeffs;
    end

    %% Getters (Boolean)
    function tf = is_constant(obj)
        % Check if polynomial is constant.
        tf = is_constant(obj.coeffs);
    end

    function tf = is_equal(obj,p)
        % Check if polynomials are equal.
        tf = is_equal@casos.package.core.GenericPolynomial(obj,p) ...
            && is_equal(obj.coeffs,p.coeffs);
    end

    function tf = is_indet(obj)
        % Check if polynomial is a vector of indeterminates.
        tf = (is_monom(obj) && is_homogeneous(obj,1));
    end

    % check if (symbolic) polynomial is linear
    tf = is_linear(obj,x);

    function tf = is_monom(obj)
        % Check if polynomial is vector of monomials (legacy).
        tf = (is_monom(obj.get_sparsity) && is_coeff_one(obj));
    end

    function tf = is_one(obj)
        % Check if polynomial is equal to identity.
        tf = (is_zerodegree(obj) && is_coeff_one(obj));
    end

    function tf = is_symbolic(obj)
        % Check if polynomial has symbolic coefficients.
        tf = is_symbolic(obj.coeffs(coeff_find(obj.get_sparsity)));
    end

    % check well-posedness
    tf = is_wellposed(obj);

    function tf = is_zero(obj)
        % Check if polynomial is equal to zero.
        tf = (is_zerodegree(obj) && is_zero(obj.coeffs));
    end

    %% Matrix & Sparsity operations
    % polynomial basis (deprecated)
    S = basis(obj,I);

    % block concatenation
    p = blockcat(a,b,c,d);

    % public RedefinesParen interface
    p = cat(dim,varargin);

    % nonzero coordinates
    [c,S] = coordinates(obj,S);

    % operator matrix (deprecated)
    [M,S] = op2basis(obj,varargin);

    % nonzero polynomial coordinates (deprecated)
    [c,S] = poly2basis(obj,varargin);

    function obj = project(obj,S)
        % Project onto sparsity pattern.
        if ~check_sz_equal(obj,S)
            throw(casos.package.core.IncompatibleSizesError.other(obj,S));
        end

        [S,obj.coeffs] = coeff_project(obj.get_sparsity,obj.coeffs,casos.Sparsity(S));
        obj = obj.set_sparsity(S);
    end

    function obj = repmat(obj,varargin)
        % Repeat polynomial matrix.
        assert(nargin > 1, 'Not enough input arguments.')

        % repition scheme
        rep = horzcat(varargin{:});

        [S,obj.coeffs] = coeff_repmat(obj.get_sparsity,obj.coeffs,rep);
        obj = obj.set_sparsity(S);
    end

    function obj = reshape(obj,varargin)
        % Reshape polynomial matrix.
        assert(length(varargin{1}) <= 2, 'Size vector must not exceed two elements.')
        assert(length(varargin) <= 2, 'Size arguments must not exceed two scalars.')

        S = reshape(obj.get_sparsity,varargin{:});
        obj = obj.set_sparsity(S);
    end

    function obj = sparsity_cast(obj,S)
        % Cast sparsity pattern onto polynomial.
        assert(nnz(obj) == nnz(S), 'Mismatching number of nonzero coefficients.')

        S = casos.Sparsity(S);

        obj.coeffs = sparsity_cast(obj.coeffs,coeff_sparsity(S));
        obj = obj.set_sparsity(S);
    end

    %% Polynomial operations
    % adjoint operator
    b = adjoint(a);

    % Frobenius norm (operators only)
    r = Fnorm2(obj)

    % polynomial integral
    b = int(a,x,varargin);
    
    % polynomial differentiation
    b = nabla(a,x);

    % squared polynomial integral norm
    r = pnorm2(obj);

    % polynomial Taylor expansion
    c = ptaylor(a,x,b,deg);

    % remove coefficients below tolerance
    b = remove_coeffs(obj,tol);

    % substitute indeterminates
    c = subs(a,x,b);

    %% Symbolic operations
    % symbolic differentiation
    b = jacobian(a,x);

    % symbolic linearization
    c = linearize(a,x,b);

    % symbolic Taylor expansion
    c = mtaylor(a,x,b,deg);

    %% Conversion
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

    % simplify coefficients
    b = simplify(obj);

    % sparsify coefficients
    b = sparsify(obj);

    %% Unary operations
    b = power(a,n);
    b = transpose(a);
    
    function p = uminus(p)
        % Unary minus.
        p.coeffs = uminus(p.coeffs);
    end

    function p = uplus(p)
        % Unary plus.
        p.coeffs = uplus(p.coeffs);
    end

    %% Binary operations
    c = compose(a,b);
    c = dot(a,b);
    c = evaluate(a,b);

    function c = minus(a,b)
        % Substract two polynomials.
        c = plus(a, uminus(b));
    end

    c = ldivide(a,p);
    c = plus(a,b);
    c = rdivide(p,b);
    c = times(a,b);

    %% Matrix operations
    c = kron(a,b);
    c = mldivide(a,p);
    b = mpower(a,n);
    c = mrdivide(p,b);
    c = mtimes(a,b);
    b = prod(a,dim);
    b = sum(a,dim);

    %% Misc
    function l = list_of_coeffs(obj)
        % Return a list of coefficients.
        l = coeff_list(obj.get_sparsity,obj.coeffs);
    end

    %% Display
    % string representation
    out = str(obj);
end

methods (Access=protected)
    % helper for static constructors
    obj = new_sym(obj,dstr,varargin);

    % protected RedefinesParen interface
    obj = parenAssign(obj,idx,varargin);
    varargout = parenReference(obj,index);
end

end
