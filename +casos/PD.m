classdef (InferiorClasses = {?casadi.DM, ?casos.Indeterminates}) ...
        PD < casos.package.core.Polynomial
% Polynomial with constant coefficients of type casadi.DM.
%
% Constructor summary:
%
%   PD()
%
% create empty (0x0) polynomial.
%
%   PD(int,int)
%
% create all-sparse polynomial.
%
%   PD(double | DM)
%
% convert double matrix.
%
%   PD(Sparsity)
%
% create from list of monomials 
% (all non-sparse coefficients equate to 1).
%
%   PX(Sparsity,scalar double | DM)
%
% create polynomial with constant coefficients 
% (all non-sparse coefficients equate to the given value).
%
%   PX(Sparsity,vector double | DM)
%
% create polynomial with constant coefficients from vector of nonzeros.
%

methods (Static,Access=protected)
    %% Polynomial interface
    function c = new_coeff(varargin)
        % New coefficients.
        c = casadi.DM(varargin{:});
    end

    function p = new_poly(varargin)
        % New polynomial.
        p = casos.PD(varargin{:});
    end
end

methods (Static)
    %% Static constructors
    function p = sym(varargin) %#ok<STOUT>
        % Symbolic variable.
        error('Symbolics not implemented for PD.');
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
    %% Conversion
    function d = casadi.DM(p)
        % Convert to double matrix.
        assert(is_zerodegree(p), 'Can only convert polynomial of degree zero.')

        d = reshape(p.coeffs,size(p));
    end

    %% Getter & symbolic framework
    function tf = is_equal(obj,p)
        % Check if polynomials are equal.
        tf = is_equal@casos.package.core.Polynomial(casos.PD(obj),casos.PD(p));
    end

    function tf = is_linear(obj,p) %#ok<STOUT,INUSD>
        % Check if polynomial is linear.
        error('Symbolics not implemented for PD.')
    end

    function varargout = jacobian(f,x) %#ok<INUSD,STOUT>
        % Symbolic Jacobian.
        if isa(x,'casos.package.core.AlgebraicObject') && is_indet(x)
            error('Use NABLA for the polynomial Jacobian.')
        end
        error('Symbolics not implemented for PD.')
    end

    function c = linearize(a,x,b) %#ok<STOUT,INUSD>
        % Symbolic linearization.
        error('Symbolics not implemented for PD.')
    end

    function c = mtaylor(a,x,b,varargin) %#ok<STOUT,INUSD>
        % Symbolic Taylor expansion.
        error('Symbolics not implemented for PD.')
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
        % Addition.
        c = plus@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function c = times(a,b)
        % Element-wise multiplication.
        c = times@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function r = dot(a,b)
        % Scalar cot product.
        r = dot@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function c = mtimes(a,b)
        % Matrix multiplication.
        c = mtimes@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function c = kron(a,b)
        % Kronecker product.
        c = kron@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function c = ldivide(a,b)
        % Element-wise left division.
        c = ldivide@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function c = rdivide(a,b)
        % Element-wise right division.
        c = rdivide@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function c = mldivide(a,b)
        % Matrix left division.
        c = mldivide@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function c = mrdivide(a,b)
        % Matrix right division.
        c = mrdivide@casos.package.core.Polynomial(casos.PD(a),casos.PD(b));
    end

    function b = cat(dim,varargin)
        % Concatenation.
        if nargin == 3
            b = cat@casos.package.core.Polynomial(dim,casos.PD(varargin{1}),casos.PD(varargin{2}));
        else
            b = cat@casos.package.core.Polynomial(dim,varargin{:});
        end
    end

    function p = blockcat(a,b,c,d)
        % Block concatenation.
        p = blockcat@casos.package.core.Polynomial(casos.PD(a),casos.PD(b),casos.PD(c),casos.PD(d));
    end

    function c = subs(a,x,b)
        % Substitution.
        c = subs@casos.package.core.Polynomial(casos.PD(a),x,casos.PD(b));
    end

    function c = ptaylor(a,x,b,deg)
        % Polynomial Taylor expansion.
        c = ptaylor@casos.package.core.Polynomial(casos.PD(a),x,casos.PD(b),deg);
    end
end

end
