classdef OperatorSparsity < casos.package.core.AbstractSparsity
% Operator sparsity class.

properties (GetAccess=private, SetAccess=immutable)
    % sparsity of linear map from in-nnz to out-nnz
    sparsity_M = casadi.Sparsity;

    % polynomial sparsity patterns of input and output
    sparsity_in = casos.Sparsity;
    sparsity_out = casos.Sparsity;
end

properties (Dependent)
    nvars;
    nterm;
    maxdeg;
    mindeg;
end

methods (Access={?casos.Sparsity, ?casos.package.core.AbstractSparsity})
    %% Friend getters
    function D = get__degsum(obj)
        % Relative degree of each monomial.
        D = reshape(obj.sparsity_out.degsum,[],1) - reshape(obj.sparsity_in.degsum,1,[]);
    end

    function get__coeffs(~)
        % Coefficient matrix.
        error('Notify the developers.')
    end

    function get__degmat(~)
        % Degree matrix.
        error('Notify the developers.')
    end

    function get__indets(~)
        % Indeterminate variables.
        error('Notify the developers.')
    end

    function get__matdim(~)
        % Matrix dimensions.
        error('Notify the developers.')
    end

    function S = get__sparsity_M(obj)
        % Sparsity of linear map.
        S = obj.sparsity_M;
    end

    function S = get__sparsity_in(obj)
        % Input sparsity pattern.
        S = obj.sparsity_in;
    end

    function S = get__sparsity_out(obj)
        % Output sparsity pattern.
        S = obj.sparsity_out;
    end
end

methods (Access=private)
    %% Private constructor
    function obj = OperatorSparsity(S,Si,So)
        % New operator sparsity pattern.
        obj.sparsity_M = casadi.Sparsity(S);
        obj.sparsity_in = Si;
        obj.sparsity_out = So;
    end
end

methods (Static)
    %% Public constructor
    function S = pattern(S,Si,So)
        % Create operator sparsity pattern.
        if nargin < 1
            % nothing to do
            error('Forbidden.')

        elseif isa(S,'casos.Sparsity')
            % convert polynomial sparsity pattern
            error('Cannot convert polynomial sparsity pattern.')

        elseif nargin < 2
            % matrix multiplication pattern
            Si = casos.Sparsity.dense(size(S,2),1);
            So = casos.Sparsity.dense(size(S,1),1);

        elseif nargin < 3
            % dual operator pattern
            Si = casos.Sparsity(Si);
            So = casos.Sparsity.dense(size(S,1),1);

        elseif nargin < 4
            % construct operator pattern
            Si = casos.Sparsity(Si);
            So = casos.Sparsity(So);

        else
            error('Undefined syntax.')
        end

        assert(size(S,2) == nnz(Si), 'Input dimensions mismatch.')
        assert(size(S,1) == nnz(So), 'Output dimensions mismatch.')

        S = casos.package.core.OperatorSparsity(S,Si,So);
    end

    function S = sparse(n,m)
        % Create all-sparse operator pattern.
        Si = casos.Sparsity(m,1);
        So = casos.Sparsity(n,1);

        % create operator with sparse coefficients
        S = casos.package.core.OperatorSparsity(casadi.Sparsity(0,0),Si,So);
    end
end

methods
    %% Copy constructor
    function S = casos.Sparsity(obj)
        % Wrap sparsity pattern.
        S = casos.Sparsity.create(obj);
    end

    %% Getters (Dependent properties)
    function d = get.maxdeg(obj)
        % Return maximum degree.
        d = max(obj.sparsity_in.maxdeg, obj.sparsity_out.maxdeg);
    end

    function d = get.mindeg(obj)
        % Return minimum degree.
        d = min(obj.sparsity_in.mindeg, obj.sparsity_out.mindeg);
    end

    function n = get.nterm(obj)
        % Return number of terms.
        n = numel(obj.sparsity_M);
    end

    function n = get.nvars(obj)
        % Return number of indeterminate variables.
        n = length(obj.indeterminates);
    end

    %% Getters
    function x = indeterminates(obj)
        % Return indeterminate variables.
        x = combine(indeterminates(obj.sparsity_in), indeterminates(obj.sparsity_out));
    end

    function n = nnz(obj)
        % Return number of nonzeros.
        n = nnz(obj.sparsity_M);
    end

    function n = numel(obj)
        % Return number of elements.
        n = numel(obj.sparsity_in)*numel(obj.sparsity_out);
    end

    function varargout = size(obj,varargin)
        % Return size of operator.
        [varargout{1:nargout}] = size(sparse(numel(obj.sparsity_out),numel(obj.sparsity_in)),varargin{:});
    end

    %% Getters (Boolean)
    function tf = is_dense(obj)
        % Check if sparsity pattern is dense.
        tf = is_dense(obj.sparsity_M);
    end

    function tf = is_dual(obj)
        % Check if operator is a linear form (dual).
        tf = (is_zerodegree(obj.sparsity_out) && isscalar(obj.sparsity_out));
    end

    function tf = is_equal(obj,op)
        % Check if operators are equal.
        tf = is_operator(op) ...
            && is_equal(obj.sparsity_in,op.sparsity_in) ...
            && is_equal(obj.sparsity_out,op.sparsity_out) ...
            && is_equal(obj.sparsity_M,op.sparsity_M);
    end

    function tf = is_matrix(obj)
        % Check if operator is mapping between vectors.
        tf = (iscolumn(obj.sparsity_in) && iscolumn(obj.sparsity_out));
    end

    function tf = is_operator(~)
        % Check for operator.
        tf = true;
    end

    function tf = is_wellposed(obj)
        % Check if operator is well formed.
        tf = is_wellposed(obj.sparsity_in) ...
            && is_wellposed(obj.sparsity_out) ...
            && size(obj.sparsity_M,1) == nnz(obj.sparsity_out) ...
            && size(obj.sparsity_M,2) == nnz(obj.sparsity_in);
    end

    function tf = is_zerodegree(obj)
        % Check if operator is of input/output degree zero.
        tf = (obj.sparsity_in.maxdeg == 0 && obj.sparsity_out.maxdeg == 0);
    end

    %% Conversion
    function S = primalize(obj)
        % Convert dual operator to primal.
        assert(is_dual(obj), 'Operator must be linear form (dual).')

        S = casos.Sparsity(obj.sparsity_in);

        assert(is_dense(obj.sparsity_M), 'Notify the developers.')
    end

    function S = dualize(obj)
        % Convert to dual.
        assert(is_dual(obj), 'Operator must be linear form (dual).')

        S = casos.Sparsity(obj);
    end

    function S = to_vector(obj,varargin)
        % Convert a scalar dual operator to a vector.
        assert(is_dual(obj) && isscalar(obj), 'Operator must be scalar linear form.')

        S = casos.Sparsity( ...
                            obj.sparsity_M, ...
                            to_vector(obj.sparsity_in,varargin{:}), ...
                            obj.sparsity_out ...
        );
    end

    function [S,I] = restrict_terms(obj,deg)
        % Restrict monomial terms.
        [Si,I1] = restrict_terms(obj.sparsity_in,deg);
        [So,I2] = restrict_terms(obj.sparsity_out,deg);

        % coefficients to keep
        I = (reshape(I2,[],1) & reshape(I1,1,[]));

        % restrict coefficients
        M = reshape(sub(obj.sparsity_M,find(I)-1),numel(So),numel(Si));

        % new sparsity pattern
        S = casos.Sparsity(M,Si,So);
    end

    %% Matrix sparsity interface
    function S = casadi.Sparsity(obj)
        % Convert zero-degree matrix pattern to casadi.Sparsity type.
        assert(is_zerodegree(obj), 'Can only convert pattern of degree zero.')
        assert(is_matrix(obj), 'Can only convert operator patterns between vectors.')

        S = casadi.Sparsity(obj.sparsity_M);
    end

    %% Display output
    function s = str(obj)
        % Return string representation.
        s = compose('[%s]->[%s],%dnz',to_char(obj.sparsity_in),to_char(obj.sparsity_out),nnz(obj.sparsity_M));
    end

    function s = signature(obj)
        % Return signature.
        in = signature(obj.sparsity_in); out = signature(obj.sparsity_out);
        s = compose('(%s)->(%s),%dnz',in{:},out{:},nnz(obj.sparsity_M));
    end

    function print_matrix(obj)
        % Print operator matrix pattern.
        disp(obj.sparsity_M)
    end
end

methods (Access={?casos.Sparsity, ?casos.package.core.AbstractSparsity})
    %% Friend class interface
    function S = coeff_sparsity(obj)
        % Return sparsity pattern of linear map.
        S = casadi.Sparsity(obj.sparsity_M);
    end

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

    % protected interface for linear operators
    [S,I1,I2] = op_intersect(obj,S2);
    [S,I1,I2] = op_join(obj,S2);
end

end
