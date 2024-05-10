classdef (InferiorClasses = {?casadi.DM, ?casos.Indeterminates}) ...
        PD < casos.package.core.Polynomial
% Polynomial with constant coefficients of type casadi.DM.

methods (Static,Access=protected)
    %% Polynomial interface
    function c = new_coeff(varargin)
        % New coefficients.
        c = casadi.DM(varargin{:});
    end

    function sym_coeff(varargin)
        % Symbolic coefficients.
        error('Not implemented.')
    end

    function p = new_poly(varargin)
        % New polynomial.
        p = casos.PD(varargin{:});
    end
end

methods (Static)
    %% Static constructors
    function sym(varargin)
        % Symbolic variable.
        error('Not implemented for PD.');
    end

    function p = empty(varargin)
        % Empty polynomial matrix.
        p = casos.PD(double.empty(varargin{:}));
    end

    function p = zeros(varargin)
        % Zero polynomial.
        p = casos.PD(casadi.DM.zeros(varargin{:}));
    end

    function p = ones(varargin)
        % One polynomial.
        p = casos.PD(casadi.DM.ones(varargin{:}));
    end

    function p = eye(varargin)
        % Identity matrix.
        p = casos.PD(casadi.DM.eye(varargin{:}));
    end
end

methods (Access=protected)
    %% Protected interface
    function tf = is_coeff_one(obj)
        % Check if nonzero coefficients are one.
        tf = is_one(obj.coeffs(coeff_find(get_sparsity(obj))));
    end
end

methods
    %% Algebraic operations
    function c = plus(a,b)
        % Addition
        c = plus@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function c = times(a,b)
        % Element-wise multiplication
        c = times@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function c = mtimes(a,b)
        % Matrix multiplication
        c = mtimes@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function b = cat(dim,varargin)
        % Concatenation.
        if nargin == 3
            b = cat@casos.package.core.Polynomial(dim,casos.PD(varargin{1}),casos.PD(varargin{2}));
        else
            b = cat@casos.package.core.Polynomial(dim,varargin{:});
        end
    end
end

end
