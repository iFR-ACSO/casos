classdef (InferiorClasses = {?casadi.Sparsity}) Sparsity ...
        < casos.package.core.PolynomialInterface
% Polynomial sparsity class.

properties (GetAccess=protected, SetAccess=private)
    % polynomial sparsity is stored in multi-index fashion, that is,
    %
    %   p = sum_a c_a*x^a
    %
    % where x = (x1,...,xN) are the indeterminates, a = (a1,...,aN) are the
    % indices/degrees, the expression x^a is shorthand for x1^a1*...*xN^aN,
    % and c_a are (possibly matrix-valued) sparsity patterns.
    coeffs = casadi.Sparsity;    % matrix of which the rows are vec(c_a)
    degmat = sparse([]);         % matrix of which the rows are [a1 ... aN]
    indets = casos.Indeterminates;             % variables {x1,...,xN} 
    matdim = [1 1];              % dimensions (size) of coefficients c_a
end

properties (Dependent=true)
    nvars;      % number of indeterminate variables
    nterm;      % number of monomial terms
    maxdeg;
    mindeg;

%     nnz;        % number of nonzero coefficients
%     matrix_nnz; % number of nonzero components
end

properties (Dependent=true, Access=protected)
    degsum;     % total degree of each monomial
end

methods
    %% Public constructor
    function obj = Sparsity(varargin)
        % Create polynomial sparsity pattern.
        if nargin == 0
            % nothing to do (null)

        elseif isa(varargin{1},'casos.Sparsity')
            % copy sparsity pattern
            assert(nargin == 1,'Too many arguments.')
            
            obj.coeffs = casadi.Sparsity(varargin{1}.coeffs);
            obj.degmat = varargin{1}.degmat;
            obj.indets = varargin{1}.indets;
            obj.matdim = varargin{1}.matdim;

        elseif ischar(varargin{1}) || isa(varargin{1},'casos.Indeterminates')
            % indeterminate (pvar / mpvar syntax)
            indets = casos.Indeterminates(varargin{:});
            N = length(indets);

            % sort variables alphabetically
            [obj.indets,~] = sort(indets);

            % return scalar linear sparsity pattern
            obj.coeffs = casadi.Sparsity.dense(N,1);
            obj.degmat = speye(N);
            obj.matdim = [1 1];

        else
            % zero-degree sparsity (casadi syntax)
            S = casadi.Sparsity(varargin{:});
            % store size and coefficients
            obj = casos.Sparsity.coeff_zerodegree(S);
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

    function d = get.degsum(obj)
        % Total degrees of each monomial.
        d = full(sum(obj.degmat,2));
    end

    function d = get.mindeg(obj)
        % Minimum degree of polynomial.
        d = min(obj.degsum);
    end

    function d = get.maxdeg(obj)
        % Maximum degree of polynomial.
        d = max(obj.degsum);
    end

    function n = nnz(obj)
        % Return number of nonzero coefficients.
        n = nnz(obj.coeffs);  
    end

    function n = numel(obj)
        % Return number of elements.
        n = prod(obj.matdim);
    end

    function n = matrix_nnz(obj)
        % Return number of nonzero (matrix) components.
        n = length(unique(get_col(obj.coeffs)));
    end

    function varargout = size(obj,varargin)
        % Return size of polynomial.
        [varargout{1:nargout}] = size(sparse(obj.matdim(1),obj.matdim(2)),varargin{:});
    end

    function x = indeterminates(obj)
        % Return indeterminate variables of polynomial.
        x = casos.Indeterminates(obj.indets);
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

    function tf = is_zerodegree(obj)
        % Check if polynomial is of degree zero.
        tf = (obj.maxdeg == 0);
    end

    function tf = is_null(obj)
        % Check if sparsity pattern is null.
        tf = is_null(obj.coeffs);
    end
end

methods (Static)
    %% Static constructors
    S = scalar(varargin);

    function S = diag(varargin)
        % Create diagonal matrix pattern.
        S = casos.Sparsity(casadi.Sparsity.diag(varargin{:}));
    end

    function S = dense(varargin)
        % Create dense matrix pattern.
        S = casos.Sparsity(casadi.Sparsity.dense(varargin{:}));
    end

    % to be completed
end

methods
    %% Polynomial Sparsity interface
    function idx = find(obj)
        % Return indices of nonzero elements.
        idx = matrix_find(obj); % TODO
    end

    function [i,j] = get_triplet(obj)
        % Return triplets for polynomial sparsity pattern.
        [i,j] = matrix_triplet(obj); % TODO
    end

    function z = monomials(obj)
        % Return scalar monomial sparsity pattern.
        z = build_monomials(obj.degmat,obj.indets);
    end

    function obj = transpose(obj), end % nothing to do
    function obj = ctranspose(obj), end % nothing to do

    function S = horzcat(varargin)
        % Horizontal concatenation.
        S = cat(2,varargin{:});
    end

    function S = vertcat(varargin)
        % Vertical concatenation.
        S = cat(1,varargin{:});
    end

    %% Conversion & matrix Sparsity interface
    function S = reshape(obj,varargin)
        % Reshape polynomial matrix.
        assert(length(varargin{1}) <= 2, 'Size vector must not exceed two elements.')
        assert(length(varargin) <= 2, 'Size arguments must not exceed two scalars.')

        S = casos.Sparsity(obj);
        S.matdim = size(reshape(sparse(size(obj,1),size(obj,2)),varargin{:}));
    end

    function S = repmat(obj,varargin)
        % Repeat sparsity pattern.
        S = coeff_repmat(obj,obj.coeffs,varargin{:});
    end

    function S = casadi.Sparsity(obj)
        % Convert zero-degree pattern to casadi.Sparsity type.
        assert(is_zerodegree(obj), 'Can only convert pattern of degree zero.')

        S = reshape(obj.coeffs,obj.matdim);
    end

    function S = matrix_sparsity(obj)
        % Return matrix sparsity pattern.
        S = reshape(sum1(obj.coeffs),obj.matdim);
    end

    function [i,j] = matrix_triplet(obj)
        % Return triplets for matrix sparsity pattern.
        [i,j] = ind2sub(size(obj),matrix_find(obj)-1);
    end

    function idx = matrix_find(obj)
        % Return indices of nonzero elements.
        idx = find(sum1(obj.coeffs));
    end 

    %% Display output
    function matrix_spy(obj)
        % Print matrix sparsity pattern.
        spy(matrix_sparsity(obj));
    end

    function disp(obj)
        % Display object.
        if is_null(obj), disp('NULL');
        else, disp@casos.package.core.PolynomialInterface(obj);
        end
    end
end

methods (Access=private)
    function obj = set_coefficients(obj,coeffs)
        % Store (generic) coefficients.
        if isa(coeffs,'casadi.Sparsity')
            obj.coeffs = coeffs;
        else
            obj.coeffs = sparsity(coeffs);
        end
    end
end

methods (Static, Access={?casos.package.core.PolynomialInterface})
    %% Static friend interface
    function [S,coeffs] = coeff_zerodegree(A)
        % Set sparsity and coefficients for zero-degree polynomial.
        S = casos.Sparsity;
        % reshape to coefficient matrix
        coeffs = reshape(A,1,numel(A));
        % store size
        S.degmat = sparse(1,0);
        S.matdim = size(A);
        % store coefficients
        S = set_coefficients(S,coeffs);
    end
end

methods (Access={?casos.package.core.PolynomialInterface})
    %% Friend class interface
    function varargout = coeff_size(obj,varargin)
        % Return size of coefficient matrix.
        [varargout{1:nargout}] = size(obj.coeffs,varargin{:});
    end

    function S = coeff_sparsity(obj)
        % Return sparsity pattern of coefficient matrix.
        S = casadi.Sparsity(obj.coeffs);
    end

    function [i,j] = coeff_triplet(obj)
        % Return triplets for coefficient matrix sparsity.
        [i,j] = get_triplet(obj.coeffs);
    end

    function idx = coeff_find(obj)
        % Return indices of nonzero coefficients.
        idx = find(obj.coeffs);
    end

    % protected interface for operations
    [S,coeffs] = coeff_repmat(obj,coeffs,varargin);
    [S,coeffs] = coeff_subsref(obj,coeffs,ii,sz);
    [S,coeffs] = coeff_subsasgn(obj,S2,coeffs,coeff2,ii);
    [S,coeffs] = coeff_transpose(obj,coeffs);
    [S,cf1,cf2] = coeff_expand(S1,S2,coeff1,coeff2);
    [S,coeffs] = coeff_cat(S1,S2,coeff1,coeff2,dim);
    [S,coeffs] = coeff_plus(S1,S2,coeff1,coeff2);
    [S,coeffs] = coeff_times(S1,S2,coeff1,coeff2);
    [S,coeffs] = coeff_mtimes(S1,S2,coeffs1,coeff2);

    % protected interface for display output
    out = str_monoms(obj,flag);

    % protected interface for subsref getters
    [monoms,L] = get_monoms(obj,I);
    [degree,L] = get_degree(obj,I);
    [indets,L] = get_indets(obj,I);

    function nv = get_nvars(obj,I)
        % Get number of variables.
        [~,L] = get_indets(obj);
        nv = nnz(any(L(I,:),1));
    end

    function nt = get_nterm(obj,I)
        % Get number of terms.
        [~,L] = get_degmat(obj);
        nt = nnz(any(L(I,:),1));
    end
end

end
