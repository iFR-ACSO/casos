classdef (InferiorClasses = {?casadi.SX, ?casadi.DM}) PS < matlab.mixin.indexing.RedefinesParen
% Polynomials with symbolic coefficients.

properties (GetAccess=protected, SetAccess=private)
    % polynomials are stored in multi-index fashion, that is,
    %
    %   p = sum_a c_a*x^a
    %
    % where x = (x1,...,xN) are the indeterminates, a = (a1,...,aN) are the
    % indices/degrees, the expression x^a is shorthand for x1^a1*...*xN^aN,
    % and c_a are (possibly matrix-valued) coefficients.
    coeffs = casadi.SX;     % matrix of which the rows are vec(c_a)
    degmat = sparse([]);    % matrix of which the rows are [a1 ... aN]
    indets = {};            % cell array of strings {x1,...,xN} 
    matdim = [0 0];         % dimensions (size) of coefficients c_a
end

properties (Dependent=true)
    nvars;      % number of indeterminate variables
    nterm;      % number of monomial terms
    maxdeg;
    mindeg;
end

methods
    %% Public constructor
    function obj = PS(varargin)
        % Create polynomial variable.
        if nargin == 0
            % nothing to do

        elseif isa(varargin{1},'casos.PS') && nargin < 2
            % keep polynomial
            obj = varargin{1};

        elseif isa(varargin{1},'casos.PS')
            % basis notation: Z'*q
            Z = varargin{1};
            q = varargin{2};

            % repeat scalar coefficient for all monomials
            if isscalar(q), q = repmat(q,size(Z,1),1); end

            assert(size(Z,1) == size(q,1), 'Incompatible size (Expected %d, got %d).', size(Z,1), size(q,1))

            % TODO: perform operation internally
            obj = Z'*q;

            if nargin > 2
                % reshape to given size
                obj.matdim = varargin{3};
            end

        elseif isa(varargin{1},'char')
            % indeterminate (pvar / mpvar syntax)
            var = varargin{1};
            arg = varargin(2:end);
            if nargin == 1
                % syntax PS('x')
                n = 1; m = 1;
                obj.indets = {var};
                iv = 1;
            elseif ischar([arg{:}])
                % syntax PS('x','y',...)
                n = length(varargin); m = 1;
                % sort variables alphabetically
                [obj.indets,iv] = unique(varargin);
            else
                % syntax PS('x',m,n)
                [n,m] = size(zeros(arg{:}));
                obj.indets = compose('%s_%d',var,1:(n*m));
                iv = 1:(n*m);
            end

            % return variables in requested order
            obj.coeffs = casadi.SX.triplet(iv-1,0:(n*m-1),ones(1,n*m),n*m,n*m);
            obj.degmat = speye(n*m);
            obj.matdim = [n m];

        else
            % constant polynomial (casadi syntax)
            C = casadi.SX(varargin{:});
            % store size and coefficients
            obj.coeffs = reshape(C,1,numel(C));
            obj.degmat = sparse(1,0);
            obj.matdim = size(C);
        end
    end

    %% Getter
    function n = get.nvars(obj)
        % Number of indeterminate variables.
        n = length(obj.indets);
    end

    function n = get.nterm(obj)
        % Number of monomials.
        n = size(obj.degmat,1);
    end

    function d = get.mindeg(obj)
        % Minimum degree of polynomial.
        d = full(min(sum(obj.degmat,2)));
    end

    function d = get.maxdeg(obj)
        % Maximum degree of polynomial.
        d = full(max(sum(obj.degmat,2)));
    end

    function varargout = size(obj,varargin)
        % Return size of polynomial.
        [varargout{1:nargout}] = size(sparse(obj.matdim(1),obj.matdim(2)),varargin{:});
    end

    function x = indeterminates(obj)
        % Return indeterminate variables of polynomial.
        x = casos.PS(obj.indets{:});
    end

    function tf = isrow(obj)
        % Check if polynomial is a row vector.
        tf = (size(obj,1) == 1);
    end

    function tf = iscolumn(obj)
        % Check if polynomial is a column vector.
        tf = (size(obj,2) == 1);
    end

    function tf = isvector(obj)
        % Check if polynomial is a vector.
        tf = any(size(obj) == 1);
    end

    function tf = isscalar(obj)
        % Check if polyonomial is a scalar.
        tf = all(size(obj) == 1);
    end

    function tf = is_symbolic(obj)
        % Check if polynomial has symbolic coefficients.
        tf = is_symbolic(obj.coeffs);
    end

    function tf = is_symexpr(obj)
        % Check if polynomial contains symbolic expressions.
        tf = ~is_constant(obj.coeffs);
    end

    function tf = is_symgram(obj)
        % Check if polynomial is in symbolic Gram form.
        tf = ~isempty(grammatrix(obj));
    end

    function tf = is_zerodegree(obj)
        % Check if polynomial is of degree zero.
        tf = (obj.maxdeg == 0);
    end

    function tf = is_zero(obj)
        % Check if polynomial is equal to zero.
        tf = (is_zerodegree(obj) && is_zero(obj.coeffs));
    end

    function tf = is_constant(obj)
        % Check if polynomial is constant.
        tf = (is_zerodegree(obj) && ~is_symexpr(obj));
    end

    function tf = is_monom(obj)
        % Check if polynomial is a vector of monomials.

        % A polynomial is a vector of monomials iff all dimensions have
        % exactly one nonzero coefficient 1
        tf = is_one(obj.coeffs == 0 | obj.coeffs == 1) && is_one(sum(obj.coeffs,1));
    end

    function tf = is_indet(obj)
        % Check if polynomial is a vector of indeterminates.

        % A polynomial is a vector of indeterminates iff it is a vector of
        % monomials with degrees equal to one.
        tf = is_monom(obj) && all(sum(obj.degmat,2) == 1);
    end
end

methods (Static)
    %% Static constructors
    p = sym(dstr,w,sz,type); % symbolic constructor

    function p = empty()
        % Create empty polynomial matrix.
        p = casos.PS;
    end

    function p = zeros(varargin)
        % Create zero polynomial.
        p = casos.PS(casadi.SX.zeros(varargin{:}));
    end

    function p = ones(varargin)
        % Create one polynomial.
        p = casos.PS(casadi.SX.ones(varargin{:}));
    end

    function p = eye(varargin)
        % Create identity matrix polynomial.
        p = casos.PS(casadi.SX.eye(varargin{:}));
    end
end

methods
    % public RedefinesParen interface
    p = cat(dim,varargin);

    function obj = reshape(obj,varargin)
        % Reshape polynomial matrix.
        assert(length(varargin{1}) <= 2, 'Size vector must not exceed two elements.')
        assert(length(varargin) <= 2, 'Size arguments must not exceed two scalars.')

        obj.matdim = size(reshape(sparse(size(obj,1),size(obj,2)),varargin{:}));
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

    %% Conversion
    function d = casadi.SX(p)
        % Convert degree-zero polynomial to casadi.SX type.
        assert(is_zerodegree(p), 'Can only convert polynomial of degree zero.')

        d = reshape(p.coeffs,p.matdim);
    end

    function d = casadi.DM(p)
        % Convert constant polynomial to casadi.DM type.
        assert(is_constant(p), 'Can only convert constant polynomial.')

        d = casadi.DM(casadi.SX(p));
    end

    function d = double(p)
        % Convert constant polynomial to double data type.
        d = full(casadi.DM(p));
    end
end

methods (Access=protected)
    % protected RedefinesParen interface
    obj = parenAssign(obj,idx,varargin);
    obj = parenDelete(obj,idx);
    n = parenListLength(obj,idx,context);
    varargout = parenReference(obj,index);

    % protected interface for subsref getters
    [monoms,L] = get_monoms(p,I);
    [degree,L] = get_degree(p,I);
    [indets,L] = get_indets(p,I);
end

end