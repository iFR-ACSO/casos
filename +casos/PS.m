classdef (InferiorClasses = {?casadi.DM, ?casadi.SX, ?casos.PD, ?casos.Indeterminates}) ...
        PS < casos.package.core.Polynomial
% Polynomial with symbolic coefficients of type casadi.SX.
%
% Constructor summary:
%
%   PS()
%
% create empty (0x0) polynomial.
%
%   PS(int,int)
%
% create all-sparse polynomial.
%
%   PS(double | DM | SX)
%
% convert double or symbolic matrix.
%
%   PS(PD)
%
% convert polynomial with constant coefficients.
%
%   PS(char, ...)
%   PS(char,int,[int])
%
% create vector (matrix) of indeterminate variables.
%
%   PS(Sparsity)
%
% create from list of monomials 
% (all non-sparse coefficients equate to 1).
%
%   PS(Sparsity,scalar double | DM)
%
% create polynomial with constant coefficients 
% (all non-sparse coefficients equate to the given value).
%
%   PS(Sparsity,scalar SX)
%
% create polynomial with symbolic coefficients
% (all non-sparse coefficients equate to the given expression).
%
%   PS(Sparsity,vector double | DM)
%
% create polynomial with constant coefficients from vector of nonzeros.
%
%   PS(Sparsity,vector SX)
%
% create polynomial with symbolic coefficients from vector of nonzeros.
%
%   PS(matrix double | DM,Sparsity,Sparsity)
%
% create operator with constant coefficients.
%
%   PS(matrix SX,Sparsity,Sparsity)
%
% create operator with symbolic coefficients.
%
% Static constructor summary:
%
%   sym(char,[vector int])
%
% create zero-degree symbolic polynomial.
%
%   sym(char,[vector int],scalar Sparsity)
%
% create symbolic polynomial (matrix) with given monomials.
%
%   sym(char,Sparsity)
%
% create symbolic polynomial or operator with given sparsity pattern.
%
%   sym(char,[...],"gram")
%
% create symbolic polynomial in Gram form.
%
%   empty([int],[int])
%   eye([int],[int])
%   ones([int],[int])
%   zeros([int],[int])
%
% create zero-degree matrix polynomial.
%
%   id_operator(int,[int])
%   one_operator(int,[int])
%   zero_operator(int,[int])
%
% create zero-degree operator.
%


methods (Static,Access=protected)
    %% Polynomial interface
    function c = new_coeff(varargin)
        % New coefficients.
        c = casadi.SX(varargin{:});
    end

    function p = new_poly(varargin)
        % New polynomial.
        p = casos.PS(varargin{:});
    end
end

methods (Static)
    %% Static constructors
    function p = sym(dstr,varargin)
        % Symbolic variable.
        S = casos.PS.sym_pattern(varargin{:});
        p = casos.PS(S,casadi.SX.sym(dstr,nnz(S)));
    end

    function p = empty(varargin)
        % Empty polynomial matrix.
        p = casos.PS(double.empty(varargin{:}));
    end

    function p = zeros(varargin)
        % Zero polynomial.
        p = casos.PS(casadi.SX.zeros(varargin{:}));
    end

    function p = ones(varargin)
        % One polynomial.
        p = casos.PS(casadi.SX.ones(varargin{:}));
    end

    function p = eye(varargin)
        % Identity matrix.
        p = casos.PS(casadi.SX.eye(varargin{:}));
    end

    function p = id_operator(varargin)
        % Identity operator.
        p = casos.PS(casos.Sparsity.diag_operator(varargin{:}),1);
    end

    function p = one_operator(varargin)
        % One operator.
        p = casos.PS(casos.Sparsity.dense_operator(varargin{:}),1);
    end

    function p = zero_operator(varargin)
        % Zero operator.
        p = casos.PS(casos.Sparsity.dense_operator(varargin{:}),0);
    end
end

methods
    %% Conversion
    function d = casadi.SX(obj)
        % Convert to double matrix.
        assert(is_zerodegree(obj), 'Can only convert polynomial of degree zero.')

        d = reshape(obj.coeffs,size(obj));
    end

    function p = casos.PD(obj)
        % Convert to double polynomial.
        assert(~is_symexpr(obj), 'Cannot convert symbolic polynomials.')

        p = casos.PD(obj.get_sparsity);
        p.coeffs = casadi.DM(obj.coeffs);
    end

    function f = to_mxfunction(obj,varargin) %#ok<STOUT,INUSD>
        % Convert to MX function.
        error('Conversion to MX not possible.')
    end

    %% Getter & symbolic framework
    function tf = is_equal(obj,p)
        % Check if polynomials are equal.
        tf = is_equal@casos.package.core.Polynomial(casos.PS(obj),casos.PS(p));
    end

    function tf = is_linear(obj,p)
        % Check if polynomial is linear.
        tf = is_linear@casos.package.core.Polynomial(casos.PS(obj),casos.PS(p));
    end

    function D = jacobian(f,x)
        % Symbolic Jacobian.
        D = jacobian@casos.package.core.Polynomial(casos.PS(f),casos.PS(x));
    end

    function c = linearize(a,x,b)
        % Symbolic linearization.
        c = linearize@casos.package.core.Polynomial(casos.PS(a),casos.PS(x),casos.PS(b));
    end

    function c = mtaylor(a,x,b,deg)
        % Symbolic Taylor expansion.
        c = mtaylor@casos.package.core.Polynomial(casos.PS(a),casos.PS(x),casos.PS(b),deg);
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
        c = plus@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function c = times(a,b)
        % Element-wise multiplication.
        c = times@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function r = dot(a,b)
        % Scalar cot product.
        r = dot@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function c = mtimes(a,b)
        % Matrix multiplication.
        c = mtimes@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function c = kron(a,b)
        % Kronecker product.
        c = kron@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function c = ldivide(a,b)
        % Element-wise left division.
        c = ldivide@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function c = rdivide(a,b)
        % Element-wise right division.
        c = rdivide@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function c = mldivide(a,b)
        % Matrix left division.
        c = mldivide@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function c = mrdivide(a,b)
        % Matrix right division.
        c = mrdivide@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function b = cat(dim,varargin)
        % Concatenation.
        if nargin == 3
            b = cat@casos.package.core.Polynomial(dim,casos.PS(varargin{1}),casos.PS(varargin{2}));
        else
            b = cat@casos.package.core.Polynomial(dim,varargin{:});
        end
    end

    function p = blockcat(a,b,c,d)
        % Block concatenation.
        p = blockcat@casos.package.core.Polynomial(casos.PS(a),casos.PS(b),casos.PS(c),casos.PS(d));
    end

    function c = subs(a,x,b)
        % Substitution.
        c = subs@casos.package.core.Polynomial(casos.PS(a),x,casos.PS(b));
    end

    function c = ptaylor(a,x,b,deg)
        % Polynomial Taylor expansion.
        c = ptaylor@casos.package.core.Polynomial(casos.PS(a),x,casos.PS(b),deg);
    end


    function c = compose(a,b)
        % Compose operators.
        c = compose@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end

    function c = evaluate(a,b)
        % Evaluate operator.
        c = evaluate@casos.package.core.Polynomial(casos.PS(a),casos.PS(b));
    end
end

end
