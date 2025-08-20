classdef (InferiorClasses = {?casadi.Sparsity, ?casadi.DM, ?casadi.SX, ?casadi.MX}) ...
    Sparsity < casos.package.core.PolynomialInterface
% Sparsity of polynomials and polynomial operators.
%
% Constructor summary:
%
%   Sparsity(int,int)
%
% create all-sparse polynomial pattern.
%
%   Sparsity(char, ...)
%   Sparsity(char,int,[int])
%
% create scalar monomial pattern of indeterminate variables.
%
%   Sparsity(matrix sparsity)
%
% create zero-degree polynomial pattern
% (equal to given matrix sparsity pattern).
%
%   Sparsity(matrix sparsity, Sparsity, Sparsity)
%
% create operator sparsity pattern.
%
% Static constructor summary:
%
%   scalar
%
% create scalar zero-degree polynomial sparsity pattern.
%
%   scalar(Indeterminates,vector int=1)
%
% create scalar monomial sparsity pattern.
%
%   diag(int,[int])
%   dense(int,[int])
%   band(int,int)
%   banded(int,int)
%   nonzeros(int,int,int)
%   triplet(int,int,int,int)
%
% create zero-degree polynomial sparsity pattern
% (equal to matrix sparsity pattern which is diagonal, dense, band, banded,
% or with nonzeros described by linear indices or triplet).
% 
%   diag(int,[int],scalar Sparsity)
%   dense(int,[int],scalar Sparsity)
%   band(int,int,scalar Sparsity)
%   banded(int,int,scalar Sparsity)
%   nonzeros(int,int,int,scalar Sparsity)
%   triplet(int,int,int,int,scalar Sparsity)
%
% create polynomial sparsity pattern
% (all non-sparse entries have given monomials).
%
%   diag(int,[int],vector Sparsity)
%   dense(int,[int],vector Sparsity)
%   band(int,int,vector Sparsity)
%   banded(int,int,vector Sparsity)
%   nonzeros(int,int,int,vector Sparsity)
%   triplet(int,int,int,int,vector Sparsity)
%
% create polynomial sparsity pattern from vector of monomial patterns.
%
%   diag_operator(int,[int])
%   dense_operator(int,[int])
%   sparse_operator(int,int)
%
% create zero-degree operator pattern.
%

properties (Access=private)
    % sparsity pattern
    pattern;
end

properties (Dependent)
    nvars;
    nterm;
    maxdeg;
    mindeg;

    sparsity_in;
    sparsity_out;
end

properties (Dependent, Access={?casos.package.core.AbstractSparsity})
    % total degree of each monomial
    degsum;
end

properties (Dependent, Access={?casos.package.core.PolynomialSparsity})
    % polynomial sparsity
    coeffs;
    degmat;
    indets;
    matdim;
end

properties (Dependent, Access={?casos.package.core.OperatorSparsity})
    % operator sparsity
    sparsity_M;
end

methods
    function obj = Sparsity(varargin)
        % Create sparsity pattern.
        if nargin == 0
            % nothing to do (null)

        elseif isa(varargin{1},'casos.Sparsity') && nargin < 2
            % copy sparsity pattern (shallow copy)
            obj.pattern = varargin{1}.pattern;

        elseif nargin == 3 && isa(varargin{3},'casos.Sparsity')
            % operator sparsity pattern
            obj.pattern = casos.package.core.OperatorSparsity.pattern(varargin{:});

        else
            % polynomial sparsity pattern
            obj.pattern = casos.package.core.PolynomialSparsity.pattern(varargin{:});
        end
    end

    %% Getters (Dependent properties)
    function C = get.coeffs(obj)
        % Coefficient matrix.
        assert(~is_null(obj), 'Null pointer.')
        C = get__coeffs(obj.pattern);
    end

    function D = get.degmat(obj)
        % Degree matrix.
        assert(~is_null(obj), 'Null pointer.')
        D = get__degmat(obj.pattern);
    end

    function d = get.degsum(obj)
        % Total degrees.
        assert(~is_null(obj), 'Null pointer.')
        d = get__degsum(obj.pattern);
    end

    function x = get.indets(obj)
        % Indeterminate variables.
        assert(~is_null(obj), 'Null pointer.')
        x = get__indets(obj.pattern);
    end

    function sz = get.matdim(obj)
        % Matrix dimensions.
        assert(~is_null(obj), 'Null pointer.')
        sz = get__matdim(obj.pattern);
    end

    function d = get.maxdeg(obj)
        % Maximum degree.
        assert(~is_null(obj), 'Null pointer.')
        d = obj.pattern.maxdeg;
    end

    function d = get.mindeg(obj)
        % Minimum degree.
        assert(~is_null(obj), 'Null pointer.')
        d = obj.pattern.mindeg;
    end

    function n = get.nterm(obj)
        % Number of monomials.
        assert(~is_null(obj), 'Null pointer.')
        n = obj.pattern.nterm;
    end

    function n = get.nvars(obj)
        % Number of indeterminate variables.
        assert(~is_null(obj), 'Null pointer.')
        n = obj.pattern.nvars;
    end

    function S = get.sparsity_in(obj)
        % Input sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        S = get__sparsity_in(obj.pattern);
    end

    function S = get.sparsity_M(obj)
        % Sparsity of linear map.
        assert(~is_null(obj), 'Null pointer.')
        S = get__sparsity_M(obj.pattern);
    end

    function S = get.sparsity_out(obj)
        % Output sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        S = get__sparsity_out(obj.pattern);
    end

    %% Getters
    function x = indeterminates(obj)
        % Indeterminate variables.
        assert(~is_null(obj), 'Null pointer.')
        x = indeterminates(obj.pattern);
    end

    function n = nnz(obj)
        % Number of nonzeros.
        assert(~is_null(obj), 'Null pointer.')
        n = nnz(obj.pattern);
    end

    function n = numel(obj)
        % Number of elements.
        assert(~is_null(obj), 'Null pointer.')
        n = numel(obj.pattern);
    end

    function varargout = size(obj,varargin)
        % Return size of pattern.
        if is_null(obj), [varargout{1:nargout}] = size(0,varargin{:});
        else, [varargout{1:nargout}] = size(obj.pattern,varargin{:});
        end
    end
    
    %% Getters (Boolean)
    function tf = is_null(obj)
        % Check if pattern is null.
        tf = isempty(obj.pattern);
    end

    function tf = is_dense(obj)
        % Check if sparsity pattern is dense.
        assert(~is_null(obj), 'Null pointer.')
        tf = is_dense(obj.pattern);
    end

    function tf = is_dual(obj)
        % Check for dual vector.
        assert(~is_null(obj), 'Null pointer.')
        tf = is_dual(obj.pattern);
    end

    function tf = is_equal(obj,S)
        % Check if sparsity patterns are equal.
        assert(~is_null(obj), 'Null pointer.')
        tf = is_equal(obj.pattern,S);
    end

    function tf = is_matrix(obj)
        % Check if operator is mapping between vectors.
        assert(~is_null(obj), 'Null pointer.')
        tf = is_matrix(obj.pattern);
    end

    function tf = is_monom(obj)
        % Check if polynomial sparsity pattern is a vector of monomials.
        assert(~is_null(obj), 'Null pointer.')
        tf = is_monom(obj.pattern);
    end

    function tf = is_operator(obj)
        % Check if pattern is operator.
        assert(~is_null(obj), 'Null pointer.')
        tf = is_operator(obj.pattern);
    end

    function tf = is_wellposed(obj)
        % Check if pattern is well-posed.
        assert(~is_null(obj), 'Null pointer.')
        tf = is_wellposed(obj.pattern);
    end

    function tf = is_zerodegree(obj)
        % Check if pattern is zero degree.
        assert(~is_null(obj), 'Null pointer.')
        tf = is_zerodegree(obj.pattern);
    end

    %% Conversion
    function S = primalize(obj)
        % Convert to primal (polynomial) pattern.
        assert(~is_null(obj), 'Null pointer.')
        S = primalize(obj.pattern);
    end

    function S = dualize(obj)
        % Convert to dual (operator) pattern.
        assert(~is_null(obj), 'Null pointer.')
        S = dualize(obj.pattern);
    end

    function S = to_vector(obj,varargin)
        % Convert to vector.
        assert(~is_null(obj), 'Null pointer.')
        S = to_vector(obj.pattern,varargin{:});
    end

    %% Concatenation
    function S = horzcat(varargin)
        % Horizontal concatenation.
        S = cat(2,varargin{:});
    end

    function S = vertcat(varargin)
        % Vertical concatenation.
        S = cat(1,varargin{:});
    end

    %% Polynomial Sparsity interface
    function idx = find(obj)
        % Return indices of nonzero elements.
        assert(~is_null(obj), 'Null pointer.')
        idx = find(obj.pattern);
    end

    function [i,j] = get_triplet(obj)
        % Return triplets for polynomial sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        [i,j] = get_triplet(obj.pattern);
    end

    function [S,varargout] = gram(obj)
        % Return Gram basis of a monomial sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        assert(~is_operator(obj), 'Not allowed for operators.')
        [S,varargout{1:nargout-1}] = gram(obj.pattern);
    end

    function [S,varargout] = grambasis(obj,varargin)
        % Return Gram basis of a polynomial sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        assert(~is_operator(obj), 'Not allowed for operators.')
        [S,varargout{1:nargout-1}] = grambasis(obj.pattern,varargin{:});
    end

    function z = monomials(obj) 
        % Return scalar monomial sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        z = monomials(obj.pattern);
    end

    function [S,I] = restrict_terms(obj,varargin)
        % Restrict terms.
        assert(~is_null(obj), 'Null pointer.')
        [S,I] = restrict_terms(obj.pattern,varargin{:});
    end

    function obj = transpose(obj), end % nothing to do
    function obj = ctranspose(obj), end % nothing to do

    %% Matrix sparsity interface
    function S = casadi.Sparsity(obj) 
        % Convert to casadi.Sparsity type.
        assert(~is_null(obj), 'Null pointer.')
        S = casadi.Sparsity(obj.pattern);
    end

    function idx = matrix_find(obj)
        % Indices of nonzero elements.
        assert(~is_null(obj), 'Null pointer.')
        idx = matrix_find(obj.pattern);
    end

    function n = matrix_nnz(obj) 
        % Return number of nonzero (matrix) components.
        assert(~is_null(obj), 'Null pointer.')
        n = matrix_nnz(obj.pattern);
    end

    function S = matrix_sparsity(obj)
        % Matrix sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        S = matrix_sparsity(obj.pattern);
    end

    function [i,j] = matrix_triplet(obj)
        % Return triplets for matrix sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        [i,j] = matrix_triplet(obj.pattern);
    end

    function S = reshape(obj,varargin) 
        % Reshape sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        S = reshape(obj.pattern,varargin{:});
    end

    %% Extended interface for CasADi data types
    function varargout = coordinates(M,S)
        % Fall back to polynomial.
        [varargout{1:nargout}] = coordinates(casos.package.polynomial(M),S);
    end

    function p = project(M,S)
        % Fall back to polynomial.
        p = project(casos.package.polynomial(M),S);
    end

    function p = sparsity_cast(M,S)
        % Fall back to polynomial.
        p = sparsity_cast(casos.package.polynomial(M),S);
    end

    function varargout = poly2basis(M,S)
        % Fall back to polynomial.
        warning('Deprecated. Use "coordinates" instead.')
        [varargout{1:nargout}] = coordinates(casos.package.polynomial(M),S);
    end

    function varargout = op2basis(M,S)
        % Fall back to operator.
        warning('Deprecated. Use "coordinates" instead.')
        [varargout{1:nargout}] = coordinates(casos.package.polynomial(M),S);
    end

    %% Display output
    function matrix_spy(obj)
        % Print matrix sparsity pattern.
        spy(matrix_sparsity(obj));
    end

    function print_terms(obj)
        % Print monomial terms.
        assert(~is_null(obj), 'Null pointer.')
        out = str_terms(obj.pattern);
        % print matrix of terms
        disp_matrix(obj,'[]',out);
        disp(' ')
    end

    function dim = signature(obj,varargin)
        % Return signature.
        assert(~is_null(obj), 'Null pointer.')
        dim = signature(obj.pattern,varargin{:});
    end

    function spy(obj)
        % Visualize polynomial sparsity pattern.
        assert(~is_null(obj), 'Null pointer.')
        assert(~is_operator(obj), 'Not allowed for operators.')
        spy(obj.pattern);
    end

    function out = str(obj)
        % String representation.
        if is_null(obj), out = {'NULL'};
        else, out = str(obj.pattern);
        end
    end

    %% Misc
    function l = list_of_degree(obj) 
        % Return a list of degrees.
        assert(~is_null(obj), 'Null pointer.')
        l = list_of_degree(obj.pattern);
    end

    function l = list_of_indets(obj) 
        % Return a list of indeterminate variables.
        assert(~is_null(obj), 'Null pointer.')
        l = list_of_indets(obj.pattern);
    end
end

methods (Static)
    %% Polynomial constructors
    function S = scalar(varargin)
        % Create scalar sparsity pattern.
        S = casos.Sparsity.create(casos.package.core.PolynomialSparsity.scalar(varargin{:}));
    end

    function S = diag(n,k,varargin)
        % Create diagonal sparsity pattern.
        if nargin < 2, S = casos.Sparsity(casadi.Sparsity.diag(n)); % diag(n)
        elseif isnumeric(k), S = casos.Sparsity(casadi.Sparsity.diag(n,k),varargin{:}); % diag(n,k,...)
        else, S = casos.Sparsity(casadi.Sparsity.diag(n),k,varargin{:}); % diag(n,...)
        end
    end

    function S = dense(n,k,varargin)
        % Create dense sparsity pattern.
        if nargin < 2, S = casos.Sparsity(casadi.Sparsity.dense(n)); % dense(n)
        elseif isnumeric(k), S = casos.Sparsity(casadi.Sparsity.dense(n,k),varargin{:}); % dense(n,k,...)
        else, S = casos.Sparsity(casadi.Sparsity.dense(n),k,varargin{:}); % dense(n,...)
        end
    end

    function S = band(n,p,varargin)
        % Create band sparsity pattern.
        S = casos.Sparsity(casadi.Sparsity.band(n,p),varargin{:});
    end

    function S = banded(n,p,varargin)
        % Create banded sparsity pattern.
        S = casos.Sparsity(casadi.Sparsity.banded(n,p),varargin{:});
    end

    function S = nonzeros(n,m,idx,varargin)
        % Create sparsity pattern with nonzeros.
        S = casos.Sparsity(casadi.Sparsity.nonzeros(n,m,idx),varargin{:});
    end

    function S = triplet(n,m,i,j,varargin)
        % Create sparsity pattern with triplets.
        S = casos.Sparsity(casadi.Sparsity.triplet(n,m,i,j),varargin{:});
    end
    
    %% Operator constructors
    function S = dense_operator(varargin)
        % Create dense operator pattern.
        S = casos.Sparsity.create(casos.package.core.OperatorSparsity.pattern(casadi.Sparsity.dense(varargin{:})));
    end

    function S = diag_operator(varargin)
        % Create diagonal operator pattern.
        S = casos.Sparsity.create(casos.package.core.OperatorSparsity.pattern(casadi.Sparsity.diag(varargin{:})));
    end

    function S = sparse_operator(varargin)
        % Create all-sparse operator pattern.
        S = casos.Sparsity.create(casos.package.core.OperatorSparsity.sparse(varargin{:}));
    end
end

methods (Static, Access={?casos.package.core.AbstractSparsity})
    %% Static friend interface
    function S = create(pattern)
        % Create new sparsity pattern.
        S = casos.Sparsity;
        S.pattern = pattern;
    end
end

methods (Access={?casos.package.core.PolynomialInterface, ?casos.package.core.AbstractSparsity})
    %% Friend class interface
    function S = coeff_sparsity(obj)
        % Return sparsity pattern of coefficients.
        assert(~is_null(obj), 'Null pointer.')
        S = coeff_sparsity(obj.pattern);
    end

    function varargout = coeff_size(obj,varargin)
        % Return size of coefficient matrix.
        assert(~is_null(obj), 'Null pointer.')
        [varargout{1:nargout}] = coeff_size(obj.pattern,varargin{:});
    end

    function [i,j] = coeff_triplet(obj)
        % Return triplets for coefficient matrix sparsity.
        assert(~is_null(obj), 'Null pointer.')
        [i,j] = coeff_triplet(obj.pattern);
    end

    function idx = coeff_find(obj)
        % Return indices of nonzero coefficients.
        assert(~is_null(obj), 'Null pointer.')
        idx = coeff_find(obj.pattern);
    end

    function cfs = coeff_repterms(obj,coeffs,nt) 
        % Repeat terms.
        assert(~is_null(obj), 'Null pointer.')
        cfs = coeff_repterms(obj.pattern,coeffs,nt);
    end

    function l = coeff_list(obj,coeffs)
        % Return list of coefficients.
        assert(~is_null(obj), 'Null pointer.')
        l = coeff_list(obj.pattern,coeffs);
    end

    %% Protected interface for polynomial operations
    function [S,coeffs] = coeff_adjoint(obj,coeffs) 
        % Coefficient matrix of adjoint.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_adjoint(obj.pattern,coeffs);
    end

    function [S,coeffs] = coeff_blkcat(obj,S2,S3,S4,cf1,cf2,cf3,cf4) 
        % Coefficients matrix of block concatenation.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_blkcat(obj.pattern,S2,S3,S4,cf1,cf2,cf3,cf4);
    end

    function [S,coeffs] = coeff_cat(obj,S2,coeff1,coeff2,dim)
        % Coefficients matrix of pairwise concatenation.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_cat(obj.pattern,S2,coeff1,coeff2,dim);
    end

    function [S,coeffs] = coeff_dot(obj,S2,coeff1,coeff2)
        % Coefficient matrix of dot product.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_dot(obj.pattern,S2,coeff1,coeff2);
    end
    
    function [cf1,cf2] = coeff_expand(obj,S2,coeff1,coeff2) 
        % Expand coefficient matrices.
        assert(~is_null(obj), 'Null pointer.')
        [cf1,cf2] = coeff_expand(obj,S2,coeff1,coeff2);
    end

    function [S,coeffs] = coeff_int(obj,coeffs,x,range) 
        % Coefficient matrix of integral.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_int(obj.pattern,coeffs,x,range);
    end

    function [S,coeffs] = coeff_kron(obj,S2,coeff1,coeff2) 
        % Coefficient matrix of Kronecker product.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_kron(obj.pattern,S2,coeff1,coeff2);
    end

    function [S,coeffs] = coeff_mtimes(obj,S2,coeffs) 
        % Coefficient matrix of matrix multiplication.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_mtimes(obj.pattern,S2,coeffs);
    end

    function [S,coeffs] = coeff_nabla(obj,coeffs,x) 
        % Coefficient matrix of partial derivative.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_nabla(obj.pattern,coeffs,x);
    end

    function [S,coeffs] = coeff_plus(obj,S2,coeff1,coeff2)
        % Coefficient matrix of addition.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_plus(obj.pattern,S2,coeff1,coeff2);
    end
    
    function [S,coeffs] = coeff_power(obj,coeffs,deg) 
        % Coefficient matrix of element-wise power.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_power(obj.pattern,coeffs,deg);
    end

    function [S,coeffs] = coeff_prod(obj,coeffs,dim) 
        % Coefficient matrix of matrix product.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_prod(obj.pattern,coeffs,dim);
    end

    function [S,coeffs] = coeff_project(obj,coeffs,S,varargin)
        % Coefficient matrix of projection.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_project(obj.pattern,coeffs,S,varargin{:});
    end

    function [S,coeffs] = coeff_repmat(obj,coeffs,varargin) 
        % Coefficient matrix of matrix repetition.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_repmat(obj.pattern,coeffs,varargin{:});
    end

    function [S,coeffs] = coeff_subsasgn(obj,S2,coeffs,coeff2,ii)
        % Coefficient matrix of subscripted assignment.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_subsasgn(obj.pattern,S2,coeffs,coeff2,ii);
    end

    function [S,coeffs] = coeff_subsref(obj,coeffs,ii,sz)
        % Coefficient matrix of subscripted reference.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_subsref(obj.pattern,coeffs,ii,sz);
    end
    
    function [S,coeffs] = coeff_substitute(obj,coeff1,x,S2,coeff2) 
        % Coefficient matrix of substitute.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_substitute(obj.pattern,coeff1,x,S2,coeff2);
    end

    function [S,coeffs] = coeff_sum(obj,coeffs,dim) 
        % Coefficient matrix of matrix sum.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_sum(obj.pattern,coeffs,dim);
    end

    function [S,coeffs] = coeff_times(obj,S2,coeff1,coeff2)
        % Coefficient matrix of element-wise product.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_times(obj.pattern,S2,coeff1,coeff2);
    end

    function [S,coeffs] = coeff_transpose(obj,coeffs) 
        % Coefficient matrix of transpose.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_transpose(obj.pattern,coeffs);
    end

    function [S,coeffs] = coeff_update(obj,coeffs,varargin) 
        % Update coefficient matrix.
        assert(~is_null(obj), 'Null pointer.')
        [S,coeffs] = coeff_update(obj.pattern,coeffs,varargin{:});
    end

    function r = coeff_properint(obj,coeffs) 
        % Coefficient matrix of proper integral.
        assert(~is_null(obj), 'Null pointer.')
        r = coeff_properint(obj.pattern,coeffs);
    end

    % Protected interface for conversion
    function x = vector_to_indeterminates(obj)
        % Convert vector of monomials to indeterminate variables.
        assert(~is_null(obj), 'Null pointer.')
        x = vector_to_indeterminates(obj.pattern);
    end

    %% Protected interface for matrix operations
    function coeffs = prepareMatrixOp(obj,coeffs,dim) 
        % Prepare for matrix operation.
        assert(~is_null(obj), 'Null pointer.')
        coeffs = prepareMatrixOp(obj.pattern,coeffs,dim);
    end

    function sz = sizeofMatrixOp(obj,dim) 
        % Return size of matrix operation.
        assert(~is_null(obj), 'Null pointer.')
        sz = sizeofMatrixOp(obj.pattern,dim);
    end

    %% Protected interface for linear operators
    function [S,varargout] = op_intersect(obj,S2)
        % Intersect polynomial sparsity patterns.
        assert(~is_null(obj), 'Null pointer.')
        [S,varargout{1:(nargout-1)}] = op_intersect(obj.pattern,S2);
    end

    function [S,varargout] = op_join(obj,S2) 
        % Join polynomial sparsity patterns.
        assert(~is_null(obj), 'Null pointer.')
        [S,varargout{1:(nargout-1)}] = op_join(obj.pattern,S2);
    end

    %% Protected interface for display output
    function out = str_monoms(obj,varargin) 
        % Return string representation for monomials.
        assert(~is_null(obj), 'Null pointer.')
        assert(~is_operator(obj), 'Not allowed for operators.')
        out = str_monoms(obj.pattern,varargin{:});
    end

    function out = str_terms(obj)
        % Return string representation for monomial terms.
        assert(~is_null(obj), 'Null pointer.')
        assert(~is_operator(obj), 'Not allowed for operators.')
        out = str_terms(obj.pattern);
    end

    % protected interface for subsref getters
    function [monoms,L] = get_monoms(obj,varargin) 
        % Return vector of monomials.
        assert(~is_null(obj), 'Null pointer.')
        [monoms,L] = get_monoms(obj.pattern,varargin{:});
    end

    function [degree,L] = get_degree(obj,varargin) 
        % Return vector of degrees.
        assert(~is_null(obj), 'Null pointer.')
        [degree,L] = get_degree(obj.pattern,varargin{:});
    end

    function [indets,L] = get_indets(obj,varargin) 
        % Return indeterminate variables.
        assert(~is_null(obj), 'Null pointer.')
        [indets,L] = get_indets(obj.pattern,varargin{:});
    end

    function nv = get_nvars(obj,varargin)
        % Get number of variables.
        assert(~is_null(obj), 'Null pointer.')
        nv = get_nvars(obj.pattern,varargin{:});
    end

    function nt = get_nterm(obj,varargin)
        % Get number of terms.
        assert(~is_null(obj), 'Null pointer.')
        nt = get_nterm(obj,varargin{:});
    end
end

methods
    %% Public interface
    function S = adjoint(obj)
        % Sparsity of adjoint operator.
        S = coeff_adjoint(obj.pattern,obj.coeff_sparsity);
    end

    function S = blockcat(S1,S2,S3,S4)
        % Block concatenation.
        S1 = casos.Sparsity(S1);
        S2 = casos.Sparsity(S2);
        S3 = casos.Sparsity(S3);
        S4 = casos.Sparsity(S4);

        % check for operators
        tf = [is_operator(S1) is_operator(S2) is_operator(S3) is_operator(S4)];

        assert(all(~tf) || all(tf), 'Must not mix polynomials and operators.')

        % block concatenate coefficient matrices
        S = coeff_blkcat(S1,S2,S3,S4,S1.coeff_sparsity,S2.coeff_sparsity,S3.coeff_sparsity,S4.coeff_sparsity);
    end

    function S = cat(dim,varargin)
        % Concatenate sparsity patterns along specified dimension.
        % Overwriting matlab.mixin.indexing.RedefinesParen.cat
        if nargin ~= 3
            % handle non-binary concatenation
            S = cat@casos.package.core.PolynomialInterface(dim,varargin{:});
            return
        end
        
        % two patterns given
        S1 = casos.Sparsity(varargin{1});
        S2 = casos.Sparsity(varargin{2});

        assert(is_operator(S1) == is_operator(S2), 'Must not mix polynomials and operators.')

        % concatenate coefficient matrices
        S = coeff_cat(S1,S2,S1.coeff_sparsity,S2.coeff_sparsity,dim);
    end

    function c = dot(a,b)
        % Dot product.
        a = casos.Sparsity(a);
        b = casos.Sparsity(b);

        if is_operator(b)
            % composition of operators
            if ~is_operator(a), a = dualize(a); end
            assert(all(size(a.sparsity_in) == size(b.sparsity_out)), 'Dimension mismatch for operator composition.')
        elseif ~is_operator(a)
            % inner product of polynomials
            assert(all(size(a) == size(b)), 'Dimension mismatch for polynomial inner product.')
        else
            % evaluation of operator on polynomial
            assert(all(size(a.sparsity_in) == size(b)), 'Dimension mismatch for operator evaluation.')
        end

        % compute dot product
        c = coeff_dot(a.pattern,b,a.coeff_sparsity,b.coeff_sparsity);
    end

    function b = integral(a,x,range)
        % Return sparsity pattern of (polynomial) integral.
        assert(~is_operator(a), 'Not allowed for operators.')
        assert(is_indet(x), 'Second argument must be vector of indeterminates.')
        
        % compute coefficient matrix of integral
        b = coeff_int(a.pattern,a.coeffs,x,range);
    end

    function c = intersect(a,b)
        % Intersect sparsity patterns.
        a = casos.Sparsity(a);
        b = casos.Sparsity(b);

        assert(is_operator(a) == is_operator(b), 'Must not mix polynomials and operators.')

        if ~check_sz_equal(a,b)
            % sparsity patterns must be of same size
            throw(casos.package.core.IncompatibleSizesError.other(a,b));
        end
        
        % intersect coefficient matrices
        c = op_intersect(a.pattern,b);
    end

    function c = kron(a,b)
        % Return sparsity pattern of Kronecker product. 
        a = casos.Sparsity(a);
        b = casos.Sparsity(b);

        assert(~is_operator(a) && ~is_operator(b), 'Not allowed for operators.')

        if isempty(a) || isempty(b)
            % handle simple case(s) for speed up:
            % product with empty polynomial is empty
            c = casos.Sparsity(size(a).*size(b));
            return
        end
        
        % else:
        % Kronecker product of coefficients
        c = coeff_kron(a.pattern,b,a.coeff_sparsity,b.coeff_sparsity);
    end

    function c = mtimes(a,b)
        % Return sparsity pattern of matrix multiplication.
        a = casos.Sparsity(a);
        b = casos.Sparsity(b);

        assert(~is_operator(a) && ~is_operator(b), 'Not allowed for operators.')

        if ~check_sz_mtimes(a,b)
            % dimensions are compatible if size(a,2) == size(b,1)
            throw(casos.package.core.IncompatibleSizesError.matrix(a,b));
        end

        % delegate to polynomial pattern
        c = mtimes(a.pattern,b);
    end

    function b = nabla(a,x)
        % Return sparsity pattern of (polynomial) Jacobian matrix.
        assert(~is_operator(a), 'Not allowed for operators.')
        assert(is_indet(x), 'Second argument must be vector of indeterminates.')
        
        % compute coefficient matrix of Jacobian
        b = coeff_nabla(a.pattern,a.coeffs,x);
    end

    function c = plus(a,b)
        % Add (join) sparsity patterns.
        a = casos.Sparsity(a);
        b = casos.Sparsity(b);

        assert(is_operator(a) == is_operator(b), 'Must not mix polynomials and operators.')

        if ~check_sz_equal(a,b)
            % sparsity patterns must be of same size
            throw(casos.package.core.IncompatibleSizesError.other(a,b));
        end
        
        % join coefficient matrices
        c = coeff_plus(a.pattern,b,a.coeff_sparsity,b.coeff_sparsity);
    end

    function b = power(a,n)
        % Return sparsity pattern of element-wise power.
        assert(~is_operator(a), 'Not allowed for operators.')

        if ~check_sz_equal(a,n)
            % sparsity pattern and exponent must be of same size
            throw(casos.package.core.IncompatibleSizesError.other(a,n));
        end
        
        % power of coefficient matrix
        b = coeff_power(a.pattern,a.coeff_sparsity,n);
    end

    function b = prod1(a)
        % Product along first dimension.
        assert(~is_operator(a), 'Not allowed for operators.')

        b = coeff_prod(a.pattern,a.coeff_sparsity,1);
    end

    function b = prod2(a)
        % Product along second dimension.
        assert(~is_operator(a), 'Not allowed for operators.')

        b = coeff_prod(a.pattern,a.coeff_sparsity,2);
    end

    function S = repmat(obj,varargin)
        % Repeat sparsity pattern.
        assert(nargin > 1, 'Not enough input arguments.')

        % repition scheme
        rep = horzcat(varargin{:});

        S = coeff_repmat(obj.pattern,obj.coeff_sparsity,rep);
    end

    function S = sub(obj,ii,S0)
        % Subreference into sparsity pattern.
        if isnumeric(S0)
            % subindices given
            [R,C] = ndgrid(ii,S0);
            ii = sub2ind(size(obj),R(:),C(:));
            sz = size(R);
        
        else
            assert(length(ii) == numel(S0), 'Not supported.')
            sz = size(S0);
        end
        
        % subreference coefficient pattern
        S = coeff_subsref(obj.pattern,obj.coeff_sparsity,ii,sz);
    end

    function c = subs(a,x,b)
        % Substitute indeterminate variables.
        a = casos.Sparsity(a);
        b = casos.Sparsity(b);

        assert(~is_operator(a),'First argument must not be an operator.')
        assert(is_indet(x),'Second argument must be vector of indeterminate variables.')
        assert(~is_operator(b),'Third argument must not be an operator.')
        assert(numel(x) == numel(b),'Second and third argument have incompatible sizes.')

        c = coeff_substitute(a.pattern,a.coeff_sparsity,x,b,b.coeff_sparsity);
    end
 
    function b = sum1(a)
        % Sum along first dimension.
        assert(~is_operator(a), 'Not allowed for operators.')

        b = coeff_sum(a.pattern,a.coeff_sparsity,1);
    end

    function b = sum2(a)
        % Sum along second dimension.
        assert(~is_operator(a), 'Not allowed for operators.')

        b = coeff_sum(a.pattern,a.coeff_sparsity,2);
    end

    function c = times(a,b)
        % Return sparsity pattern of element-wise multiplication.
        a = casos.Sparsity(a);
        b = casos.Sparsity(b);

        assert(~is_operator(a) && ~is_operator(b), 'Not allowed for operators.')

        if ~check_sz_equal(a,b)
            % sparsity patterns must be of same size
            throw(casos.package.core.IncompatibleSizesError.other(a,b));
        end
        
        % join coefficient matrices
        c = coeff_times(a.pattern,b,a.coeff_sparsity,b.coeff_sparsity);
    end

    function b = T(a) 
        % Return sparsity pattern of transpose.
        assert(~is_operator(a), 'Not allowed for operators. Use "adjoint" instead.')

        b = coeff_transpose(a,a.coeff_sparsity);
    end
end

end
