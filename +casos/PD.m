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

    function p = empty()
        % Empty polynomial matrix.
        p = casos.PD;
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

end
