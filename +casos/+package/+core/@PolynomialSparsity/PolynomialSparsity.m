classdef PolynomialSparsity < casos.package.core.AbstractSparsity
% Polynomial sparsity class.

properties (GetAccess=private, SetAccess=immutable)
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

properties (Dependent)
    nvars;      % number of indeterminate variables
    nterm;      % number of monomial terms
    maxdeg;
    mindeg;

%     nnz;        % number of nonzero coefficients
%     matrix_nnz; % number of nonzero components
end

properties (Dependent, Access=private)
    degsum;     % total degree of each monomial
end

methods (Access={?casos.Sparsity, ?casos.package.core.AbstractSparsity})
    %% Friend getters
    function D = get__degsum(obj)
        % Total degree of each monomial
        D = obj.degsum;
    end

    function C = get__coeffs(obj)
        % Coefficient matrix.
        C = obj.coeffs;
    end

    function D = get__degmat(obj)
        % Degree matrix.
        D = obj.degmat;
    end

    function x = get__indets(obj)
        % Indeterminate variables.
        x = obj.indets;
    end

    function sz = get__matdim(obj)
        % Matrix dimensions.
        sz = obj.matdim;
    end

    function get__sparsity_M(~)
        % Sparsity of linear map.
        error('Notify the developers.')
    end

    function get__sparsity_in(~)
        % Input sparsity pattern.
        error('Notify the developers.')
    end

    function get__sparsity_out(~)
        % Output sparsity pattern.
        error('Notify the developers.')
    end
end

methods (Access=private)
    %% Private constructor
    function obj = PolynomialSparsity(coeffs,degmat,indets,matdim)
        % New polynomial sparsity pattern.
        obj.coeffs = coeffs;
        obj.degmat = degmat;
        obj.indets = indets;
        obj.matdim = matdim;
    end
end

methods (Static)
    %% Public constructor
    S = scalar(varargin);

    function S = pattern(varargin)
        % Create polynomial sparsity pattern.
        if nargin == 0
            % nothing to do (null)

        elseif isa(varargin{1},'casos.Sparsity') && nargin < 2
            % copy sparsity pattern
            coeffs = casadi.Sparsity(varargin{1}.coeffs);
            degmat = varargin{1}.degmat;
            indets = varargin{1}.indets;
            matdim = varargin{1}.matdim;

        elseif ischar(varargin{1}) || isa(varargin{1},'casos.Indeterminates')
            % indeterminate (pvar / mpvar syntax)
            indets = casos.Indeterminates(varargin{:});
            N = length(indets);

            % sort variables alphabetically
            [indets,~] = sort(indets);

            % return scalar linear sparsity pattern
            coeffs = casadi.Sparsity.dense(N,1);
            degmat = speye(N);
            matdim = [1 1];

        elseif nargin == 2 && isa(varargin{2},'casos.Sparsity')
            % apply sparsity pattern to (vector of) monomials
            S = casadi.Sparsity(varargin{1});
            w = varargin{2};

            assert(~is_operator(w), 'Not allowed for operators.')

            if isscalar(w)
                % repeat monomials
                coeffs = repmat(reshape(S,1,S.numel),w.nterm,1);
            else
                % assign monomials to nonzeros
                assert(nnz(S) == numel(w),'Dimension mismatch.')

                idx = find(S);
                [ii,jj] = coeff_triplet(w);
                coeffs = casadi.Sparsity.triplet(w.nterm,S.numel,ii,idx(jj+1)-1);
            end

            degmat = w.degmat;
            indets = w.indets;
            matdim = size(S);

        else
            % zero-degree sparsity (casadi syntax)
            S = casadi.Sparsity(varargin{:});
            % reshape to coefficient matrix
            coeffs = reshape(S,1,numel(S));
            degmat = sparse(1,0);
            indets = casos.Indeterminates;
            matdim = size(S);
        end

        S = casos.package.core.PolynomialSparsity(coeffs,degmat,indets,matdim);
    end
end

methods
    %% Copy constructor
    function S = casos.Sparsity(obj,varargin)
        % Wrap sparsity pattern.
        S = casos.Sparsity.create(obj);
    end

    %% Getters (Dependent properties)
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

    %% Getters
    function x = indeterminates(obj)
        % Return indeterminate variables of polynomial.
        x = casos.Indeterminates(obj.indets);
    end

    function n = matrix_nnz(obj)
        % Return number of nonzero (matrix) components.
        n = length(unique(get_col(obj.coeffs)));
    end

    function n = nnz(obj)
        % Return number of nonzero coefficients.
        n = nnz(obj.coeffs);  
    end

    function n = numel(obj)
        % Return number of elements.
        n = prod(obj.matdim);
    end

    function varargout = size(obj,varargin)
        % Return size of polynomial.
        [varargout{1:nargout}] = size(sparse(obj.matdim(1),obj.matdim(2)),varargin{:});
    end

    %% Getters (Boolean)
    function tf = is_dense(obj)
        % Check if sparsity pattern has no sparse coefficients.
        tf = is_dense(obj.coeffs);
    end

    function tf = is_equal(obj,S)
        % Check if sparsity patterns are equal.
        tf = ~is_operator(S) ...
            && is_equal(obj.indets,S.indets) ...
            && isequal(obj.degmat,S.degmat) ...
            && is_equal(obj.coeffs,S.coeffs) ...
            && isequal(obj.matdim,S.matdim);
    end

    % check if vector monomials
    [tf,I] = is_monom(obj);

    function tf = is_operator(~)
        % Check for operator.
        tf = false;
    end

    % check well-posedness
    tf = is_wellposed(obj);

    function tf = is_zerodegree(obj)
        % Check if polynomial is of degree zero.
        tf = (obj.maxdeg == 0);
    end

    %% Conversion
    function S = primalize(obj)
        % Convert to primal.
        S = casos.Sparsity.create(obj);
    end

    function S = dualize(obj)
        % Convert to dual operator.
        S = casos.Sparsity(casadi.Sparsity.dense(1,nnz(obj)),casos.Sparsity(obj),casos.Sparsity.scalar);
    end

    S = to_vector(obj,varargin);

    %% Polynomial Sparsity interface
    % gram basis from monomials
    varargout = gram(obj);

    % gram basis for polynomial
    varargout = grambasis(obj,varargin);

    % gram unit basis
    Z = gramunit(obj);

    function z = monomials(obj)
        % Return scalar monomial sparsity pattern.
        z = build_monomials(obj.degmat,obj.indets);
    end

    % pattern of matrix multiplication
    S = mtimes(a,b);

    % restrict terms to degrees
    [S,I] = restrict_terms(obj,deg);

    %% Conversion & matrix Sparsity interface
    function idx = matrix_find(obj)
        % Return indices of nonzero elements.
        idx = find(sum1(obj.coeffs));
    end 

    function S = matrix_sparsity(obj)
        % Return matrix sparsity pattern.
        S = reshape(sum1(obj.coeffs),obj.matdim);
    end

    function [i,j] = matrix_triplet(obj)
        % Return triplets for matrix sparsity pattern.
        [i,j] = ind2sub(size(obj),matrix_find(obj)-1);
    end

    function S = reshape(obj,varargin)
        % Reshape polynomial matrix.
        assert(length(varargin{1}) <= 2, 'Size vector must not exceed two elements.')
        assert(length(varargin) <= 2, 'Size arguments must not exceed two scalars.')

        sz = size(reshape(sparse(size(obj,1),size(obj,2)),varargin{:}));
        S = casos.Sparsity(casos.package.core.PolynomialSparsity(obj.coeffs,obj.degmat,obj.indets,sz));
    end

    function S = casadi.Sparsity(obj)
        % Convert zero-degree pattern to casadi.Sparsity type.
        assert(is_zerodegree(obj), 'Can only convert pattern of degree zero.')

        S = reshape(obj.coeffs,obj.matdim);
    end

    % visualize sparsity pattern
    spy(obj);

    %% Misc
    function l = list_of_degree(obj)
        % Return a list of degrees.
        l = arrayfun(@(i) full(obj.degmat(i,:)),1:obj.nterm,'UniformOutput',false);
    end

    function l = list_of_indets(obj)
        % Return a list of indeterminate variables.
        l = str(obj.indets);
    end

    %% Display
    % signature representation
    dim = signature(obj,flag);

    function s = size_to_char(obj)
        % Return string representation of size.
        out = str(obj);
        assert(isscalar(out), 'Notify the developers.')
        s = out{1};
    end

    % string representation
    out = str(obj);
end

methods (Access={?casos.Sparsity, ?casos.package.core.AbstractSparsity})
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

    function cfs = coeff_repterms(obj,coeffs,nt)
        % Repeat terms.
        cfs = repmat(coeffs,nt/obj.nterm,1);
    end

    function l = coeff_list(obj,coeffs)
        % Return list of coefficients.
        l = arrayfun(@(i) reshape(coeffs(i,:),size(obj)),1:obj.nterm,'UniformOutput',false);
    end

    % protected interface for polynomial operations
    [S,coeffs] = coeff_adjoint(obj,coeffs);
    [S,coeffs] = coeff_blkcat(obj,S2,S3,S4,cf1,cf2,cf3,cf4);
    [S,coeffs] = coeff_cat(obj,S2,coeff1,coeff2,dim);
    [S,coeffs] = coeff_dot(obj,S2,coeff1,coeff2);
    [cf1,cf2]  = coeff_expand(obj,S2,coeff1,coeff2);
    [S,coeffs] = coeff_int(obj,coeffs,x,range);
    [S,coeffs] = coeff_kron(obj,S2,coeff1,coeff2);
    [S,coeffs] = coeff_mtimes(obj,S2,coeffs);
    [S,coeffs] = coeff_nabla(obj,coeffs,x);
    [S,coeffs] = coeff_plus(obj,S2,coeff1,coeff2);
    [S,coeffs] = coeff_power(obj,coeffs,deg);
    [S,coeffs] = coeff_prod(obj,coeffs,dim);
    [S,coeffs] = coeff_project(obj,coeffs,S,keep_zeros);
            r  = coeff_properint(obj,coeffs);
    [S,coeffs] = coeff_repmat(obj,coeffs,varargin);
    [S,coeffs] = coeff_subsref(obj,coeffs,ii,sz);
    [S,coeffs] = coeff_subsasgn(obj,S2,coeffs,coeff2,ii);
    [S,coeffs] = coeff_substitute(obj,coeff1,x,S2,coeff2)
    [S,coeffs] = coeff_sum(obj,coeffs,dim);
    [S,coeffs] = coeff_times(obj,S2,coeff1,coeff2);
    [S,coeffs] = coeff_transpose(obj,coeffs);
    [S,coeffs] = coeff_update(obj,coeffs,sz,dim);

    % manipulate nonzero coefficients
    v = coeff_getnz(obj,I);
    [S,coeffs] = coeff_setnz(obj,I,v);

    % protected interface for conversion
    x = vector_to_indeterminates(obj);

    % protected interface for matrix operations
    coeffs = prepareMatrixOp(obj,coeffs,dim);
    sz = sizeofMatrixOp(obj,dim);

    % protected interface for linear operators
    [S,I1,I2] = op_intersect(obj,S2);
    [S,I1,I2] = op_join(obj,S2);
    [S,I] = op_subsref(obj,ii);

    % protected interface for display output
    out = str_monoms(obj,flag);
    out = str_terms(obj);

    % protected interface for subsref getters
    [monoms,L] = get_monoms(obj,I);
    [degree,L] = get_degree(obj,I);
    [indets,L] = get_indets(obj,I);
end

end
