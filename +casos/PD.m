classdef (InferiorClasses = {?casadi.DM}) PD < casos.package.core.Polynomial
% Polynomial with constant coefficients of type casadi.DM

methods (Static,Access=protected)
    % implement polynomial interface
    function c = new_coeff(varargin)
        % New coefficients.
        c = casadi.DM(varargin{:});
    end

    function s = sym_coeff(varargin)
        % Symbolic coefficients.
        error('Not implemented.')
    end

    function p = new_poly(varargin)
        % New polynomial.
        p = casos.PD(varargin{:});
    end
end

end
