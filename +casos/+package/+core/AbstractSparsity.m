classdef AbstractSparsity < casos.package.core.PolynomialInterface
% Abstract sparsity class.

properties (Abstract)
    nvars;
    nterm;
    maxdeg;
    mindeg;
end

properties (Abstract, Access=protected)
    degsum;
end

methods (Abstract)
    %% Getter
    tf = is_null(obj);
    tf = is_polynomial(obj);
    tf = is_operator(obj);
    tf = is_wellposed(obj);
    tf = is_zerodegree(obj);

    n = nnz(obj);
    n = numel(obj);
    
    x = indeterminates(obj);

    %% Conversion
    S = primalize(obj);
    S = dualize(obj);

    S = to_vector(obj,I,row_vector);

    [S,I] = restrict_terms(obj,deg);

    %% Matrix sparsity interface
    S = matrix_sparsity(obj);
    idx = matrix_find(obj);

    %% Display output
    dim = signature(obj);
end

methods
    %% Common interface
    function tf = is_primal(obj)
        % Check for primal vector.
        tf = is_polynomial(obj);
    end

    function tf = is_dual(obj) %#ok<MANU>
        % Check for dual vector.
        tf = false;
    end

    function tf = is_matrix(obj) %#ok<STOUT>
        % Check if operator is mapping between vectors.
        error('Function "is_matrix" not supported by class "%s".', class(obj))
    end

    function tf = isempty(obj)
        % Check if pattern is empty.
        tf = (numel(obj) == 0);
    end

    function tf = isrow(obj)
        % Check if operator is a row vector.
        tf = (size(obj,1) == 1);
    end

    function tf = iscolumn(obj)
        % Check if operator is a column vector.
        tf = (size(obj,2) == 1);
    end

    function tf = isvector(obj)
        % Check if operator is a vector.
        tf = any(size(obj) == 1);
    end

    function tf = isscalar(obj)
        % Check if operator is a scalar.
        tf = all(size(obj) == 1);
    end

    function tf = is_equal(obj,S)
        % Check if sparsity patterns are equal.
        tf = isa(obj,class(S));
    end

    %% Polynomial Sparsity interface
    function idx = find(obj)
        % Return indices of nonzero elements.
        idx = matrix_find(obj); % TODO
    end

    function [i,j] = get_triplet(obj)
        % Return triplets for polynomial sparsity pattern.
        [i,j] = matrix_triplet(obj); % TODO
    end

    function z = monomials(obj) %#ok<STOUT>
        % Return scalar monomial sparsity pattern.
        error('Function "monomials" not supported by class "%s".', class(obj))
    end

    function S = horzcat(varargin)
        % Horizontal concatenation.
        S = cat(2,varargin{:});
    end

    function S = vertcat(varargin)
        % Vertical concatenation.
        S = cat(1,varargin{:});
    end

    function S = blockcat(S1,S2,S3,S4)
        % Block concatenation.
        S = coeff_blkcat(S1,S2,S3,S4,S1.coeff_sparsity,S2.coeff_sparsity,S3.coeff_sparsity,S4.coeff_sparsity);
    end

    function S = sum1(S)
        % Sum along first dimension.
        S = coeff_sum(S,S.coeff_sparsity,1);
    end

    function S = sum2(S)
        % Sum along second dimension.
        S = coeff_sum(S,S.coeff_sparsity,2);
    end

    function S = prod1(S)
        % Product along first dimension.
        S = coeff_prod(S,S.coeff_sparsity,1);
    end

    function S = prod2(S)
        % Product along second dimension.
        S = coeff_prod(S,S.coeff_sparsity,2);
    end

    function S = subs(S1,x,S2)
        % Substitute indeterminate variables.
        assert(is_indet(x),'Second argument must be vector of indeterminate variables.')
        assert(numel(x) == numel(S2),'Second and third argument have incompatible sizes.')

        S = coeff_subs(S1,S1.coeff_sparsity,x,S2,S2.coeff_sparsity);
    end
 
    function obj = transpose(obj), end % nothing to do
    function obj = ctranspose(obj), end % nothing to do

    function S = T(obj) 
        % Return sparsity pattern of transpose.
        S = coeff_transpose(obj,obj.coeff_sparsity);
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

    function [i,j] = matrix_triplet(obj)
        % Return triplets for matrix sparsity pattern.
        [i,j] = ind2sub(size(obj),matrix_find(obj)-1);
    end

    function n = matrix_nnz(obj) %#ok<STOUT>
        % Return number of nonzero (matrix) components.
        error('Function "matrix_nnz" not supported by class "%s".', class(obj))
    end

    %% Operator interface
    function S = adjoint(obj)
        % Sparsity of adjoint.
        if is_primal(obj)
            % adjoint is dual
            S = dualize(obj);

        elseif is_dual(obj)
            % adjoint is primal
            S = primalize(obj);

        else
            % compute adjoint
            S = coeff_adjoint(obj,obj.coeff_sparsity);
        end
    end

    %% Extended interface
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
        [varargout{1:nargout}] = poly2basis(casos.package.polynomial(M),S);
    end

    function varargout = op2basis(M,S)
        % Fall back to operator.
        [varargout{1:nargout}] = op2basis(casos.package.polynomial(M),S);
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

methods (Abstract, Access={?casos.package.core.PolynomialInterface})
    %% Friend class interface
    S = new_sparse(varargin);
    S = coeff_sparsity(obj);
    
    % protected interface for polynomial operations
    [S,coeffs] = coeff_project(obj,coeffs,S,keep_zeros);
    [S,coeffs] = coeff_subsref(obj,coeffs,ii,sz);
    [S,coeffs] = coeff_subsasgn(obj,S2,coeffs,coeff2,ii);
    [S,coeffs] = coeff_cat(S1,S2,coeff1,coeff2,dim);
    [S,coeffs] = coeff_plus(S1,S2,coeff1,coeff2);
    [S,coeffs] = coeff_times(S1,S2,coeff1,coeff2);
end

methods (Access={?casos.package.core.PolynomialInterface})
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
    function [S,coeffs] = coeff_adjoint(obj,coeffs) %#ok<STOUT,INUSD>
        % Coefficient matrix of adjoint.
        error('Function "coeff_adjoint" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_blkcat(obj,S2,S3,S4,cf1,cf2,cf3,cf4) %#ok<STOUT,INUSD>
        % Coefficients matrix of block concatenation.
        error('Function "coeff_blkcat" not supported by class "%s".', class(obj))
    end

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

    function [S,coeffs] = coeff_repmat(obj,coeffs,varargin) %#ok<STOUT,INUSD>
        % Coefficient matrix of matrix repetition.
        error('Function "coeff_repmat" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_sum(obj,coeffs,dim) %#ok<STOUT,INUSD>
        % Coefficient matrix of matrix sum.
        error('Function "coeff_sum" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_transpose(obj,coeffs) %#ok<STOUT,INUSD>
        % Coefficient matrix of transpose.
        error('Function "coeff_transpose" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_substitute(obj,coeff1,x,S2,coeff2) %#ok<STOUT,INUSD>
        % Coefficient matrix of substitute.
        error('Function "coeff_substitute" not supported by class "%s".', class(obj))
    end

    function [S,coeffs] = coeff_update(obj,coeffs,sz,dim) %#ok<STOUT,INUSD>
        % Update coefficient matrix.
        error('Function "coeff_update" not supported by class "%s".', class(obj))
    end

    function r = coeff_properint(obj,coeffs) %#ok<STOUT,INUSD>
        % Coefficient matrix of proper integral.
        error('Function "coeff_properint" not supported by class "%s".', class(obj))
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

methods
    %% Public interface
    function S = cat(dim,varargin)
        % Concatenate sparsity patterns along specified dimension.
        % Overwriting matlab.mixin.indexing.RedefinesParen.cat
        if nargin ~= 3
            % handle non-binary concatenation
            S = cat@casos.package.core.PolynomialInterface(dim,varargin{:});
            return
        end
        
        % two patterns given
        S1 = varargin{1};
        S2 = varargin{2};

        % concatenate coefficient matrices
        S = coeff_cat(S1,S2,S1.coeff_sparsity,S2.coeff_sparsity,dim);

    end

    function c = intersect(a,b)
        % Intersect sparsity patterns.
        if ~check_sz_equal(a,b)
            % sparsity patterns must be of same size
            throw(casos.package.core.IncompatibleSizesError.other(a,b));
        end
        
        % intersect coefficient matrices
        c = op_intersect(a,b);
    end

    function c = kron(a,b)
        % Return sparsity pattern of Kronecker product.        
        if isempty(a) || isempty(b)
            % handle simple case(s) for speed up:
            % product with empty polynomial is empty
            c = a.new_sparse(size(a).*size(b));
            return
        end
        
        % else:
        % Kronecker product of coefficients
        c = coeff_kron(a,b,a.coeff_sparsity,b.coeff_sparsity);
    end


    function c = plus(a,b)
        % Add (join) sparsity patterns.
        if ~check_sz_equal(a,b)
            % sparsity patterns must be of same size
            throw(casos.package.core.IncompatibleSizesError.other(a,b));
        end
        
        % join coefficient matrices
        c = coeff_plus(a,b,a.coeff_sparsity,b.coeff_sparsity);
    end

    function b = power(a,n)
        % Return sparsity pattern of element-wise power.
        if ~check_sz_equal(a,n)
            % sparsity pattern and exponent must be of same size
            throw(casos.package.core.IncompatibleSizesError.other(a,n));
        end
        
        % power of coefficient matrix
        b = coeff_power(a,a.coeff_sparsity,n);
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
        S = coeff_subsref(obj,obj.coeff_sparsity,ii,sz);
    end

    function c = times(a,b)
        % Return sparsity pattern of element-wise multiplication.
        if ~check_sz_equal(a,b)
            % sparsity patterns must be of same size
            throw(casos.package.core.IncompatibleSizesError.other(a,b));
        end
        
        % join coefficient matrices
        c = coeff_times(a,b,a.coeff_sparsity,b.coeff_sparsity);
    end
end

end
