classdef (Abstract) AbstractOperator < casos.package.core.PolynomialInterface
% A linear operator between polynomials.

properties (Access=protected)
    % linear map from in-nnz to out-nnz
    matrix;
end

properties (SetAccess=private,GetAccess=public)
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
    function obj = AbstractOperator(M,Si,So)
        % New operator.
        if nargin < 1
            % empty operator
            Si = casos.Sparsity(0,0);
            So = casos.Sparsity(0,0);
            M = [];

        elseif isa(M,'casos.package.core.AbstractOperator')
            % copy or convert operator
            Si = casos.Sparsity(M.sparsity_in);
            So = casos.Sparsity(M.sparsity_out);
            M = M.matrix;
        % Note: Conversion from polynomial is ambiguous
        % elseif isa(M,'casos.package.core.GenericPolynomial')
        %     % scalar multiplication with polynomial
        %     assert(nargin < 2,'Too many input arguments.')
        %     Zi = casos.Sparsity(1,1);
        %     [M,Zo] = poly2basis(M);

        elseif nargin < 2
            % matrix multiplication
            Si = casos.Sparsity.dense(size(M,2),1);
            So = casos.Sparsity.dense(size(M,1),1);

        elseif nargin < 3
            % dual operator
            Si = casos.Sparsity(Si);
            So = casos.Sparsity.dense(size(M,1),1);

        elseif nargin < 4
            % construct operator
            Si = casos.Sparsity(Si);
            So = casos.Sparsity(So);

        else
            error('Undefined syntax.')
        end

        assert(size(M,2) == nnz(Si), 'Input dimensions mismatch.')
        assert(size(M,1) == nnz(So), 'Output dimensions mismatch.')

        obj.matrix = obj.new_matrix(M);
        obj.sparsity_in = Si;
        obj.sparsity_out = So;
    end

    %% Getter
    function varargout = size(obj,varargin)
        % Return size of operator.
        [varargout{1:nargout}] = size(sparse(obj.numel_out,obj.numel_in),varargin{:});
    end

    function sz = size_in(obj,varargin)
        % Return input size.
        sz = size(obj.sparsity_in,varargin{:});
    end

    function sz = size_out(obj,varargin)
        % Return output size.
        sz = size(obj.sparsity_out,varargin{:});
    end

    function n = numel_in(obj)
        % Return number of input elements.
        n = numel(obj.sparsity_in);
    end

    function n = numel_out(obj)
        % Return number of output elements.
        n = numel(obj.sparsity_out);
    end

    function n = nnz_in(obj)
        % Return number input nonzeros.
        n = nnz(obj.sparsity_in);
    end

    function n = nnz_out(obj)
        % Return number of output nonzeros.
        n = nnz(obj.sparsity_out);
    end

    function n = nterm_in(obj)
        % Return number of input terms.
        n = nterm(obj.sparsity_in);
    end

    function n = nterm_out(obj)
        % Return number of output terms.
        n = nterm_out(obj.sparsity_out);
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

    function tf = is_symbolic(obj)
        % Check if operator has symbolic coefficients.
        tf = is_symbolic(obj.matrix);
    end

    function tf = is_zerodegree(obj)
        % Check if operator is of input/output degree zero.
        tf = (obj.sparsity_in.maxdeg == 0 && obj.sparsity_out.maxdeg == 0);
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
        tf = (iscolumn(obj.sparsity_in) && iscolumn(obj.sparsity_out));
    end

    function tf = is_polynomial(obj)
        % Check if operator corresponds to multiplication with polynomial.
        tf = (is_zerodegree(obj.sparsity_in) && obj.numel_in == 1);
    end

    function tf = is_dual(obj)
        % Check if operator is a linear form (dual).
        tf = (is_zerodegree(obj.sparsity_out) && obj.numel_out == 1);
    end

    function tf = is_equal(obj,op)
        % Check if operators are equal.
        tf = is_equal(obj.sparsity_in,op.sparsity_in) ...
            && is_equal(obj.sparsity_out,op.sparsity_out) ...
            && is_equal(obj.matrix,op.matrix);
    end

    function tf = is_wellposed(obj)
        % Check if operator is well formed.
        tf = is_wellformed(obj.sparsity_in) ...
            && is_wellformed(obj.sparsity_out) ...
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
