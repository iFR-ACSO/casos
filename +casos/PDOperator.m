classdef PDOperator < casos.package.core.AbstractOperator
% A linear operator between polynomials with constant coefficients.

methods (Static,Access=protected)
    %% Operator interface
    function op = new_operator(varargin)
        % New operator.
        op = casos.PDOperator(varargin{:});
    end

    function mx = new_matrix(varargin)
        % New matrix.
        mx = casadi.DM(varargin{:});
    end
end

methods (Static)
    %% Static constructors
    function op = empty(varargin)
        % Empty operator.
        op = casos.PDOperator(double.empty(varargin{:}));
    end

    function op = zeros(varargin)
        % Zero operator.
        op = casos.PDOperator(casadi.DM.zeros(varargin{:}));
    end

    function op = ones(varargin)
        % One operator.
        op = casos.PDOperator(casadi.DM.ones(varargin{:}));
    end

    function op = eye(varargin)
        % Identity operator.
        op = casos.PDOperator(casadi.DM.eye(varargin{:}));
    end
end

methods
    %% Conversion
    function D = casadi.DM(obj)
        % Convert to double matrix.
        assert(is_matrix(obj) && is_zerodegree(obj), 'Can only convert matrix operatos of degree zero.')

        I1 = find(obj.sparsity_in); I2 = find(obj.sparsity_out);
        [ii,jj] = get_triplet(sparsity(obj.matrix));
        Sp = casadi.Sparsity.triplet(size(obj,1),size(obj,2),I2(ii),I1(jj));
        D = sparsity_cast(a.matrix,Sp);
    end
    
    %% Operator algebra
    function c = dot(a,b)
        % Composition or evaluation.
        c = dot@casos.package.core.AbstractOperator(casos.PDOperator(a),b);
    end

    function c = plus(a,b)
        % Addition.
        c = plus@casos.package.core.AbstractOperator(casos.PDOperator(a),casos.PDOperator(b));
    end

    function b = cat(dim,varargin)
        % Concatenation.
        if nargin == 3
            b = cat@casos.package.core.AbstractOperator(dim,casos.PDOperator(varargin{1}),casos.PDOperator(varargin{2}));
        else
            b = cat@casos.package.core.AbstractOperator(dim,varargin{:});
        end
    end
end

end
