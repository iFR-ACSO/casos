classdef (Abstract, InferiorClasses = {?casadi.DM, ?casadi.SX, ?casadi.MX}) ...
        AlgebraicObject < matlab.mixin.indexing.RedefinesParen
% Base class for algebraic objects.

methods (Abstract)
    % check if object represents indeterminate variables
    tf = is_indet(obj);
end

methods
    %% Algebraic interface
    function t = cat(N,varargin)
        % Concatenate algebraic objects.
        args = cellfun(@(v) casos.package.polynomial(v), varargin, 'UniformOutput',false);
        t = cat(N,args{:});
    end

    function b = ctranspose(a)
        % Conjugate transpose.
        b = transpose(a);
    end

    function c = mtimes(a,b)
        % Multiply algebraic objects (matrix product).
        c = mtimes(casos.package.polynomial(a),casos.package.polynomial(b));
    end

    function c = plus(a,b)
        % Add algebraic objects.
        c = plus(casos.package.polynomial(a),casos.package.polynomial(b));
    end

    function c = times(a,b)
        % Multiply algebraic objects (element-wise).
        c = times(casos.package.polynomial(a),casos.package.polynomial(b));
    end

    function b = uminus(a)
        % Invert algebraic object.
        b = uminus(casos.package.polynomial(a));
    end

    function b = uplus(a)
        % Unary plus.
        b = uplus(casos.package.polynomial(a));
    end

    %% Common indeterminate variables interface
    function z = monomials(obj,varargin)
        % Return monomial sparsity pattern.
        assert(is_indet(obj), 'First input must be indeterminate variables.')

        z = casos.Sparsity.scalar(obj,varargin{:});
    end
end

methods (Access=protected)
    %% protected RedefinesParen interface
    function n = parenListLength(obj,indexOp,context)
        % Return length of indexing operation into algebraic object.
        if length(indexOp) < 2, n = 1; else, n = listLength(obj,indexOp(2:end),context); end
    end

    function varargout = parenDelete(~,varargin) %#ok<STOUT>
        % Overwriting matlab.mixin.indexing.RedefinesParen.parenDelete
        error('Not supported.')
    end
end

end
