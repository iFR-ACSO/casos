classdef (InferiorClasses = {?casadi.Sparsity, ?casadi.DM, ?casadi.SX, ?casadi.MX}) ...
    OperatorSparsity < casos.package.core.PolynomialInterface
% Operator sparsity class.

properties (Access=protected)
    % sparsity of linear map from in-nnz to out-nnz
    sparsity_M = casadi.Sparsity;
end

properties (SetAccess=private,GetAccess=public)
    % polynomial sparsity patterns of input and output
    sparsity_in = casos.Sparsity;
    sparsity_out = casos.Sparsity;
end

methods
    %% Public constructor
    function obj = OperatorSparsity(S,Si,So)
        % New operator sparsity pattern.
        if nargin < 1
            % nothing to do (null)

        elseif isa(S,'casos.package.core.OperatorSparsity')
            % copy operator pattern
            Si = casos.Sparsity(S.sparsity_in);
            So = casos.Sparsity(S.sparsity_out);
            S = S.sparsity_M;

        elseif nargin < 2
            % matrix multiplication pattern
            Si = casos.Sparsity.dense(size(M,2),1);
            So = casos.Sparsity.dense(size(M,1),1);

        elseif nargin < 3
            % dual operator pattern
            Si = casos.Sparsity(Si);
            So = casos.Sparsity.dense(size(M,1),1);

        elseif nargin < 4
            % construct operator pattern
            Si = casos.Sparsity(Si);
            So = casos.Sparsity(So);

        else
            error('Undefined syntax.')
        end

        assert(size(S,2) == nnz(Si), 'Input dimensions mismatch.')
        assert(size(S,1) == nnz(So), 'Output dimensions mismatch.')

        obj.sparsity_M = casadi.Sparsity(S);
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
        n = nterm(obj.sparsity_out);
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

    function tf = is_zerodegree(obj)
        % Check if operator is of input/output degree zero.
        tf = (obj.sparsity_in.maxdeg == 0 && obj.sparsity_out.maxdeg == 0);
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
            && is_equal(obj.sparsity_M,op.sparsity_M);
    end

    function tf = is_wellposed(obj)
        % Check if operator is well formed.
        tf = is_wellformed(obj.sparsity_in) ...
            && is_wellformed(obj.sparsity_out) ...
            && size(obj.sparsity_M,1) == obj.numel_out ...
            && size(obj.sparsity_M,2) == obj.numel_in;
    end

    %% Display output
    function s = str(obj)
        % Return string representation.
        s = compose('[%s]->[%s],%dnz',to_char(obj.sparsity_in),to_char(obj.sparsity_out),nnz(obj.sparsity_M));
    end

    function print_matrix(obj)
        % Print operator matrix pattern.
        disp(obj.sparsity_M)
    end
end

methods (Access={?casos.package.core.PolynomialInterface})
    %% Friend class interface
    function S = matrix_sparsity(obj)
        % Return sparsity pattern of linear map.
        S = casadi.Sparsity(obj.sparsity_M);
    end
end

end
