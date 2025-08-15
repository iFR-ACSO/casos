classdef AbstractSparsity < handle
% Abstract sparsity class.

properties (Abstract)
    nvars;
    nterm;
    maxdeg;
    mindeg;
end

methods (Abstract)
    %% Getters
    tf = is_dense(obj,S);
    tf = is_equal(obj,S);
    tf = is_operator(obj);
    tf = is_wellposed(obj);
    tf = is_zerodegree(obj);

    x = indeterminates(obj);

    n = nnz(obj);
    n = numel(obj);
    
    %% Conversion
    S = primalize(obj);
    S = dualize(obj);

    S = to_vector(obj,I,row_vector);

    [S,I] = restrict_terms(obj,deg);

    %% Display output
    dim = signature(obj);
end

methods (Abstract, Access={?casos.Sparsity, ?casos.package.core.AbstractSparsity})
    %% Friend getters
    D = get__degsum(obj);

    % polynomial sparsity
    C = get__coeffs(obj);
    D = get__degmat(obj);
    x = get__indets(obj);
    sz = get__matdim(obj);

    % operator sparsity
    S = get__sparsity_M(obj);
    S = get__sparsity_in(obj);
    S = get__sparsity_out(obj);
end

methods
    %% Common interface
    function tf = is_dual(obj) %#ok<MANU>
        % Check for dual vector.
        tf = false;
    end

    function tf = is_matrix(obj) %#ok<STOUT>
        % Check if operator is mapping between vectors.
        error('Function "is_matrix" not supported by class "%s".', class(obj))
    end

    function tf = is_monom(obj) %#ok<STOUT>
        % Check if polynomial sparsity pattern is a vector of monomials.
        error('Function "is_monom" not supported by class "%s".', class(obj))
    end

    %% Polynomial Sparsity interface
    function idx = find(obj)
        % Return indices of nonzero elements.
        idx = matrix_find(obj); % TODO
    end

    function [S,varargout] = gram(obj) %#ok<STOUT>
        % Return gram basis for a monomial sparsity pattern.
        error('Function "gram" not supported by class "%s".', class(obj));
    end

    function [S,varargout] = grambasis(obj,varargin) %#ok<STOUT>
        % Return gram basis for a polynomial sparsity pattern.
        error('Function "grambasis" not supported by class "%s".', class(obj));
    end

    function [i,j] = get_triplet(obj)
        % Return triplets for polynomial sparsity pattern.
        [i,j] = matrix_triplet(obj); % TODO
    end

    function z = monomials(obj) %#ok<STOUT>
        % Return scalar monomial sparsity pattern.
        error('Function "monomials" not supported by class "%s".', class(obj))
    end

    %% Matrix sparsity interface
    function S = casadi.Sparsity(obj) %#ok<STOUT>
        % Convert to casadi.Sparsity type.
        error('Conversion to casadi.Sparsity not supported by class "%s".', class(obj))
    end

    function S = reshape(obj,varargin) %#ok<STOUT>
        % Reshape sparsity pattern.
        error('Function "reshape" not supported by class "%s".', class(obj))
    end

    function S = matrix_sparsity(obj) %#ok<STOUT>
        % Return matrix sparsity pattern.
        error('Function "matrix_sparsity" not supported by class "%s".', class(obj))
    end

    function idx = matrix_find(obj) %#ok<STOUT>
        % Return indices of nonzero matrix elements.
        error('Function "matrix_find" not supported by class "%s".', class(obj))
    end

    function [i,j] = matrix_triplet(obj) %#ok<STOUT>
        % Return triplets for matrix sparsity pattern.
        error('Function "matrix_triplet" not supported by class "%s".', class(obj))
    end

    function n = matrix_nnz(obj) %#ok<STOUT>
        % Return number of nonzero (matrix) components.
        error('Function "matrix_nnz" not supported by class "%s".', class(obj))
    end

    %% Display output
    function print_terms(obj)
        % Print monomial terms.
        error('Function "print_terms" not supported by class "%s".', class(obj))
    end

    function spy(obj)
        % Visualize polynomial sparsity pattern.
        error('Function "spy" not supported by class "%s".', class(obj))
    end

    %% Misc
    function l = list_of_degree(obj) %#ok<STOUT>
        % Return a list of degrees.
        error('Function "list_of_degree" not supported by class "%s".', class(obj))
    end

    function l = list_of_indets(obj) %#ok<STOUT>
        % Return a list of indeterminate variables.
        error('Function "list_of_indets" not supported by class "%s".', class(obj))
    end
end

methods (Abstract, Access={?casos.Sparsity, ?casos.package.core.AbstractSparsity})
    %% Friend class interface
    S = coeff_sparsity(obj);
    
    % protected interface for polynomial operations
    [S,coeffs] = coeff_adjoint(obj,coeffs);
    [S,coeffs] = coeff_blkcat(obj,S2,S3,S4,cf1,cf2,cf3,cf4);
    [S,coeffs] = coeff_cat(obj,S2,coeff1,coeff2,dim);
    [S,coeffs] = coeff_dot(obj,S2,coeff1,coeff2);
    [S,coeffs] = coeff_plus(obj,S2,coeff1,coeff2);
    [S,coeffs] = coeff_project(obj,coeffs,S,keep_zeros);
    [S,coeffs] = coeff_repmat(obj,coeffs,varargin);
    [S,coeffs] = coeff_subsref(obj,coeffs,ii,sz);
    [S,coeffs] = coeff_subsasgn(obj,S2,coeffs,coeff2,ii);
    [S,coeffs] = coeff_update(obj,coeffs,varargin);
end

methods (Access={?casos.Sparsity, ?casos.package.core.AbstractSparsity})
    %% Friend class interface
    function varargout = coeff_size(obj,varargin)
        % Return size of coefficient matrix.
        [varargout{1:nargout}] = size(obj.coeff_sparsity,varargin{:});
    end

    function [i,j] = coeff_triplet(obj)
        % Return triplets for coefficient matrix sparsity.
        [i,j] = get_triplet(obj.coeff_sparsity);
    end

    function idx = coeff_find(obj)
        % Return indices of nonzero coefficients.
        idx = find(obj.coeff_sparsity);
    end

    function cfs = coeff_repterms(obj,coeffs,nt) %#ok<STOUT,INUSD>
        % Repeat terms.
        error('Function "coeff_repterms" not supported by class "%s".', class(obj))
    end

    function l = coeff_list(obj,coeffs) %#ok<STOUT,INUSD>
        % Return list of coefficients.
        error('Function "coeff_list" not supported by class "%s".', class(obj))
    end

    %% Protected interface for polynomial operations
    function [cf1,cf2] = coeff_expand(obj,S2,coeff1,coeff2) %#ok<STOUT,INUSD>
        % Expanded coefficient matrices.
        error('Function "coeff_expand" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_int(obj,coeffs,x,range) %#ok<STOUT,INUSD>
        % Coefficient matrix of integral.
        error('Function "coeff_int" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_kron(obj,S2,coeff1,coeff2) %#ok<STOUT,INUSD>
        % Coefficient matrix of Kronecker product.
        error('Function "coeff_kron" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_mtimes(obj,S2,coeffs) %#ok<STOUT,INUSD>
        % Coefficient matrix of matrix multiplication.
        error('Function "coeff_mtimes" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_nabla(obj,coeffs,x) %#ok<STOUT,INUSD>
        % Coefficient matrix of partial derivative.
        error('Function "coeff_nabla" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_power(obj,coeffs,deg) %#ok<STOUT,INUSD>
        % Coefficient matrix of element-wise power.
        error('Function "coeff_power" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_prod(obj,coeffs,dim) %#ok<STOUT,INUSD>
        % Coefficient matrix of matrix product.
        error('Function "coeff_prod" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_sum(obj,coeffs,dim) %#ok<STOUT,INUSD>
        % Coefficient matrix of matrix sum.
        error('Function "coeff_sum" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_times(obj,S2,coeff1,coeff2) %#ok<STOUT,INUSD>
        % Coefficient matrix of element-wise multiplication.
        error('Function "coeff_times" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_transpose(obj,coeffs) %#ok<STOUT,INUSD>
        % Coefficient matrix of transpose.
        error('Function "coeff_transpose" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_substitute(obj,coeff1,x,S2,coeff2) %#ok<STOUT,INUSD>
        % Coefficient matrix of substitute.
        error('Function "coeff_substitute" not supported by class "%s".', class(obj))
    end

    function r = coeff_properint(obj,coeffs) %#ok<STOUT,INUSD>
        % Coefficient matrix of proper integral.
        error('Function "coeff_properint" not supported by class "%s".', class(obj))
    end

    %% Protected interface for conversion
    function x = vector_to_indeterminates(obj) %#ok<STOUT>
        % Convert vector of monomials to indeterminate variables.
        error('Function "vector_to_indeterminates" not supported by class "%s".', class(obj))
    end

    %% Protected interface for matrix operations
    function coeffs = prepareMatrixOp(obj,coeffs,dim) %#ok<INUSD>
        % Prepare for matrix operation.
        error('Function "prepareMatrixOp" not supported by class "%s".', class(obj))
    end

    function sz = sizeofMatrixOp(obj,dim) %#ok<STOUT,INUSD>
        % Return size of matrix operation.
        error('Function "sizeofMatrixOp" not supported by class "%s".', class(obj))
    end

    %% Protected interface for linear operators
    function [S,I1,I2] = op_intersect(obj,S2) %#ok<STOUT,INUSD>
        % Intersect polynomial sparsity patterns.
        error('Function "op_intersect" not supported by class "%s".', class(obj))
    end

    function [S,I1,I2] = op_join(obj,S2) %#ok<STOUT,INUSD>
        % Join polynomial sparsity patterns.
        error('Function "op_join" not supported by class "%s".', class(obj))
    end

    %% Protected interface for display output
    function out = str_monoms(obj,flag) %#ok<STOUT,INUSD>
        % Return string representation for monomials.
        error('Function "str_monoms" not supported by class "%s".', class(obj))
    end

    function out = str_terms(obj) %#ok<STOUT>
        % Return string representation for terms.
        error('Function "str_terms" not supported by class "%s".', class(obj))
    end

    % protected interface for subsref getters
    function [monoms,L] = get_monoms(obj,I) %#ok<STOUT,INUSD>
        % Return vector of monomials.
        error('Function "get_monoms" not supported by class "%s".', class(obj))
    end

    function [degree,L] = get_degree(obj,I) %#ok<STOUT,INUSD>
        % Return vector of degrees.
        error('Function "get_degree" not supported by class "%s".', class(obj))
    end

    function [indets,L] = get_indets(obj,I) %#ok<STOUT,INUSD>
        % Return indeterminate variables.
        error('Function "get_indets" not supported by class "%s".', class(obj))
    end

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
