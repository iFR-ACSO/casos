classdef PSFunction < casos.package.functions.FunctionInternal
% Polynomial (PS) function interface.

properties (Access=private)
    func;

    sparsity_i;
    sparsity_o;
end

properties (SetAccess=private)
    class_name = 'PSFunction';
end

properties (Constant,Access=protected)
    allow_eval_on_basis = true;
end

methods
    function obj = PSFunction(name, ex_i, ex_o, name_i, name_o, varargin)
        % Create new casadi function object.
        obj@casos.package.functions.FunctionInternal(name);

        % parse polynomial expressions
        [sym_i,obj.sparsity_i] = cellfun(@parse_expr, ex_i, 'UniformOutput', false);
        [sym_o,obj.sparsity_o] = cellfun(@parse_expr, ex_o, 'UniformOutput', false);

        % define function between coefficients
        obj.func = casadi.Function(name, sym_i, sym_o, name_i, name_o, varargin{:});
    end

    %% Implement FunctionInternal
    function n = get_n_in(obj)
        % Number of inputs.
        n = n_in(obj.func);
    end

    function str = get_name_in(obj,varargin)
        % Name of inputs.
        str = name_in(obj.func,varargin{:});
    end

    function z = get_sparsity_in(obj,idx)
        % Monomials of inputs.
        z = obj.sparsity_i{idx+1};
    end

    function idx = get_index_in(obj,str)
        % Index of inputs.
        idx = index_in(obj.func,str);
    end

    function n = get_n_out(obj)
        % Number of outputs.
        n = n_out(obj.func);
    end

    function str = get_name_out(obj,varargin)
        % Name of outputs.
        str = name_out(obj.func,varargin{:});
    end

    function z = get_sparsity_out(obj,idx)
        % Monomials of outputs.
        z = obj.sparsity_o{idx+1};
    end

    function idx = get_index_out(obj,str)
        % Index of outputs.
        idx = index_out(obj.func,str);
    end
end

methods (Access=protected)
    %% Internal evaluation
    function argout = eval_on_basis(obj, argin)
        % Evaluate function on nonzero coordinates.
        argout = call(obj.func, argin);
    end
end

end

function [coeffs,z] = parse_expr(p)
% Return coefficients, monomials, and size of polynomial expression.

    p = casos.PS(p);
    [coeffs,z] = poly2basis(p);
end
