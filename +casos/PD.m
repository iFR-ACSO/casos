classdef (InferiorClasses = {?casadi.DM}) PD < casos.package.core.Polynomial
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
end

end
