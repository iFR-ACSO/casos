classdef PSOperator < casos.package.core.AbstractOperator
% A linear operator between polynomials with symbolic coefficients.

methods (Static,Access=protected)
    %% Operator interface
    function op = new_operator(varargin)
        % New operator.
        op = casos.PSOperator(varargin{:});
    end

    function mx = new_matrix(varargin)
        % New matrix.
        mx = casadi.SX(varargin{:});
    end
end

methods (Static)
    %% Static constructors
    function op = sym(dstr,Si,So)
        % Symbolic operator.
        if nargin < 2
            % symbolic scalar multiplication
            op = casos.PSOperator(casadi.SX.sym(dstr));

        elseif nargin < 3
            % symbolic dual operator
            op = casos.PSOperator(casadi.SX.sym(dstr,1,nnz(Si)),Si);

        else
            % symbolic operator
            op = casos.PSOperator(casadi.SX.sym(dstr,nnz(So),nnz(Si)),Si,So);
        end
    end

    function op = empty(varargin)
        % Empty operator.
        op = casos.PSOperator(double.empty(varargin{:}));
    end

    function op = zeros(varargin)
        % Zero operator.
        op = casos.PSOperator(casadi.SX.zeros(varargin{:}));
    end

    function op = ones(varargin)
        % One operator.
        op = casos.PSOperator(casadi.SX.ones(varargin{:}));
    end

    function op = eye(varargin)
        % Identity operator.
        op = casos.PSOperator(casadi.SX.eye(varargin{:}));
    end
end

methods
    %% Conversion
    function D = casadi.SX(obj)
        % Convert to double matrix.
        assert(is_matrix(obj) && is_zerodegree(obj), 'Can only convert matrix operatos of degree zero.')

        I1 = find(obj.sparsity_in); I2 = find(obj.sparsity_out);
        [ii,jj] = get_triplet(sparsity(obj.matrix));
        Sp = casadi.Sparsity.triplet(size(obj,1),size(obj,2),I2(ii),I1(jj));
        D = sparsity_cast(a.matrix,Sp);
    end

    function op = casos.PDOperator(obj)
        % Convert to double operator.
        assert(~is_symexpr(obj), 'Cannot convert symbolic operators.')

        op = casos.PDOperator(obj.matrix,obj.sparsity_in,obj.sparsity_out);
    end
    
    %% Operator algebra
    function c = dot(a,b)
        % Composition or evaluation.
        c = dot@casos.package.core.AbstractOperator(casos.PSOperator(a),b);
    end

    function c = plus(a,b)
        % Addition.
        c = plus@casos.package.core.AbstractOperator(casos.PSOperator(a),casos.PSOperator(b));
    end

    function b = cat(dim,varargin)
        % Concatenation.
        if nargin == 3
            b = cat@casos.package.core.AbstractOperator(dim,casos.PSOperator(varargin{1}),casos.PSOperator(varargin{2}));
        else
            b = cat@casos.package.core.AbstractOperator(dim,varargin{:});
        end
    end
end

end
