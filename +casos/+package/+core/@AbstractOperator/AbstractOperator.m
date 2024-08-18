classdef (Abstract) AbstractOperator < casos.package.core.PolynomialInterface
% A linear operator between polynomials.

properties (Access=protected)
    % linear map from in-nnz to out-nnz
    matrix;

    op_sparsity;    % operator sparsity pattern
end

properties (Dependent,GetAccess=public)
    % polynomial sparsity patterns of input and output
    sparsity_in;
    sparsity_out;
end

methods (Static,Abstract,Access=protected)
    op = new_operator(varargin);
    mx = new_matrix(varargin);
end

methods (Static,Abstract)
    op = empty(varargin);
    op = zeros(varargin);
    op = ones(varargin);
    op = eye(varargin);
end

methods
    %% Public constructor
    function obj = AbstractOperator(M,varargin)
        % New operator.
        if nargin < 1
            % empty operator
            S = casos.package.core.OperatorSparsity(casadi.Sparsity(0,0));
            M = [];

        elseif isa(M,'casos.package.core.AbstractOperator')
            % copy or convert operator
            S = casos.package.core.OperatorSparsity(M.op_sparsity);
            M = M.matrix;
        
        elseif isa(M,'casos.package.core.OperatorSparsity')
            % sparsity pattern syntax
            assert(nargin < 3, 'Too many arguments.')

            S = M;
            M = obj.new_matrix(matrix_sparsity(S),varargin{:});

        else
            S = casos.package.core.OperatorSparsity(sparsity(M),varargin{:});
        end

        assert(size(M,2) == nnz_in(S), 'Input dimensions mismatch.')
        assert(size(M,1) == nnz_out(S), 'Output dimensions mismatch.')

        obj.matrix = obj.new_matrix(M);
        obj.op_sparsity = S;
    end

    %% Getter
    function Si = get.sparsity_in(obj)
        % Return input polynomial sparsity pattern.
        Si = obj.op_sparsity.sparsity_in;
    end

    function So = get.sparsity_out(obj)
        % Return output polynomial sparsity pattern.
        So = obj.op_sparsity.sparsity_out;
    end

    function varargout = size(obj,varargin)
        % Return size of operator.
        [varargout{1:nargout}] = size(obj.op_sparsity,varargin{:});
    end

    function sz = size_in(obj,varargin)
        % Return input size.
        sz = size_in(obj.op_sparsity,varargin{:});
    end

    function sz = size_out(obj,varargin)
        % Return output size.
        sz = size_out(obj.op_sparsity,varargin{:});
    end

    function n = numel_in(obj)
        % Return number of input elements.
        n = numel_in(obj.op_sparsity);
    end

    function n = numel_out(obj)
        % Return number of output elements.
        n = numel_out(obj.op_sparsity);
    end

    function n = nnz_in(obj)
        % Return number input nonzeros.
        n = nnz_in(obj.op_sparsity);
    end

    function n = nnz_out(obj)
        % Return number of output nonzeros.
        n = nnz_out(obj.op_sparsity);
    end

    function n = nterm_in(obj)
        % Return number of input terms.
        n = nterm_in(obj.op_sparsity);
    end

    function n = nterm_out(obj)
        % Return number of output terms.
        n = nterm_out(obj.op_sparsity);
    end

    function S = sparsity(obj)
        % Return (copy of) operator sparsity pattern.
        S = casos.package.core.OperatorSparsity(obj.op_sparsity);
    end

    function tf = isrow(obj)
        % Check if operator is a row vector.
        tf = isrow(obj.op_sparsity);
    end

    function tf = iscolumn(obj)
        % Check if operator is a column vector.
        tf = iscolumn(obj.op_sparsity);
    end

    function tf = isvector(obj)
        % Check if operator is a vector.
        tf = isvector(obj.op_sparsity);
    end

    function tf = isscalar(obj)
        % Check if operator is a scalar.
        tf = isscalar(obj.op_sparsity);
    end

    function tf = is_symbolic(obj)
        % Check if operator has symbolic coefficients.
        tf = is_symbolic(obj.matrix);
    end

    function tf = is_zerodegree(obj)
        % Check if operator is of input/output degree zero.
        tf = is_zerodegree(obj.op_sparsity);
    end

    function tf = is_symexpr(obj)
        % Check if operator contains symbolic expressions.
        tf = ~is_constant(obj.matrix);
    end

    function tf = is_zero(obj)
        % Check if operator is equal to zero.
        tf = (is_zerodegree(obj) && is_zero(obj.matrix));
    end

    function tf = is_constant(obj)
        % Check if operator is constant.
        tf = (is_zerodegree(obj) && ~is_symexpr(obj));
    end

    function tf = is_matrix(obj)
        % Check if operator is mapping between vectors.
        tf = is_matrix(obj.op_sparsity);
    end

    function tf = is_polynomial(obj)
        % Check if operator corresponds to multiplication with polynomial.
        tf = is_polynomial(obj.op_sparsity);
    end

    function tf = is_dual(obj)
        % Check if operator is a linear form (dual).
        tf = is_dual(obj.op_sparsity);
    end

    function tf = is_equal(obj,op)
        % Check if operators are equal.
        tf = is_equal(obj.op_sparsity,op.op_sparsity) ...
            && is_equal(obj.matrix,op.matrix);
    end

    function tf = is_wellposed(obj)
        % Check if operator is well formed.
        tf = is_wellformed(obj.op_sparsity) ...
            && size(obj.matrix,1) == obj.numel_out ...
            && size(obj.matrix,2) == obj.numel_in;
    end
end

methods (Static)
    %% Static conversion
    function op = from_primal(p)
        % Convert from polynomial p such that dot(op,x) = dot(p,x).
        [P,Z] = poly2basis(p);
        op = casos.package.operator(P',Z);
    end
end

methods
    %% Algebraic operations
    function c = mtimes(a,b)
        % Multiplication with scalar.
        if isa(a,'casos.package.core.AbstractOperator')
            % right-hand scalar multiplication
            assert(isscalar(b),'Only scalar multiplication allowed.')
            c = a.new_operator(b*a.matrix,a.sparsity_in,a.sparsity_out);

        else
            % left-hand scalar multiplication
            assert(isscalar(a),'Only scalar multiplication allowed.')
            c = a.new_operator(a*b.matrix,b.sparsity_in,b.sparsity_out);
        end
    end

    function b = adjoint(a)
        % Compute adjoint operator.
        b = a.new_operator(a.matrix',a.sparsity_out,a.sparsity_in);
    end

    function b = uminus(a)
        % Negate operator.
        b = a.new_operator(uminus(a.matrix),a.sparsity_in,a.sparsity_out);
    end

    function c = minus(a,b)
        % Substract two operators.
        c = plus(a,uminus(b));
    end

    %% Conversion
    function p = to_polynomial(obj)
        % Convert to polynomial p such that p = dot(op,1).
        assert(is_polynomial(obj), 'Conversion to polynomial not possible.')
        p = casos.package.polynomial(obj.sparsity_out,obj.matrix');
    end

    function p = to_primal(obj)
        % Convert to polynomial p such that dot(op,x) = dot(p,x).
        assert(is_dual(obj), 'Conversion to primal not possible.')
        p = casos.package.polynomial(obj.sparsity_in,obj.matrix);
    end

    %% Display output
    function print_matrix(obj)
        % Print operator matrix.
        disp(obj.matrix)
    end
end

end
