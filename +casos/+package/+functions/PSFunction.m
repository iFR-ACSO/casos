classdef PSFunction < casos.package.functions.FunctionInterface
% Polynomial (PS) functionn interface.

properties (Access=private)
    func;

    monom_i;
    monom_o;

    size_i;
    size_o;
end

properties (SetAccess=protected)
    class_name = 'PSFunction';
end

methods
    function obj = PSFunction(name, ex_i, ex_o, name_i, name_o, varargin)
        % Create new casadi function object.
        obj@casos.package.functions.FunctionInterface(name);

        % parse polynomial expressions
        [sym_i,obj.monom_i,obj.size_i] = cellfun(@parse_expr, ex_i, 'UniformOutput', false);
        [sym_o,obj.monom_o,obj.size_o] = cellfun(@parse_expr, ex_o, 'UniformOutput', false);

        % define function between coefficients
        obj.func = casadi.Function(name, sym_i, sym_o, name_i, name_o, varargin{:});
    end

    %% Implement FunctionInterface
    function argout = call(obj, argin)
        % Evaluate casadi function object.
        in = cellfun(@(p,z) poly2basis(casos.PS(p),z), argin, obj.monom_i, 'UniformOutput', false);

        % call casadi function
        out = call(obj.func, in);

        % return result
        argout = cellfun(@(c,z,sz) reshape(c'*z,sz), out, obj.monom_o, obj.size_o, 'UniformOutput', false);
    end

    function n = get_n_in(obj)
        % Number of inputs.
        n = n_in(obj.func);
    end

    function str = get_name_in(obj,varargin)
        % Name of inputs.
        str = name_in(obj.func,varargin{:});
    end

    function z = get_monomials_in(obj,idx)
        % Monomials of inputs.
        z = obj.monom_i{idx+1};
    end

    function val = get_default_in(obj,idx)
        % Default inputs.
        z = get_monomials_in(idx);
        sz = get_size_in(idx);
        % build default polynomial
        val = reshape(default_in(obj.func,idx)*z,sz);
    end

    function sz = get_size_in(obj,idx)
        % Size of inputs.
        sz = obj.size_i{idx+1};
    end

    function n = get_n_out(obj)
        % Number of outputs.
        n = n_out(obj.func);
    end

    function str = get_name_out(obj,varargin)
        % Name of outputs.
        str = name_out(obj.func,varargin{:});
    end

    function z = get_monomials_out(obj,idx)
        % Monomials of outputs.
        z = obj.monom_o{idx+1};
    end

    function sz = get_size_out(obj,idx)
        % Size of outputs.
        sz = obj.size_o{idx+1};
    end
end

end

function [coeffs,z,sz] = parse_expr(p)
% Return coefficients, monomials, and size of polynomial expression.

    p = casos.PS(p);
    [coeffs,z] = poly2basis(p);
    sz = size(p);
end
