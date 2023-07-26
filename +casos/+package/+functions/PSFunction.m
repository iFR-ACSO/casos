classdef PSFunction < casos.package.functions.FunctionInterface
% Polynomial (PS) functionn interface.

properties (Access=private)
    func;
end

properties (SetAccess=protected)
    class_name = 'PSFunction';
end

properties (Dependent)
    name;
end

methods
    function obj = PSFunction(name, args_i, args_o, varargin)
        % Create new casadi function object.
        [sym_i,n_i] = cellfun(@(arg) deal(get_symbolic(arg),arg.name), args_i, 'UniformOutput', false);
        [sym_o,n_o] = cellfun(@(arg) deal(get_symbolic(arg),arg.name), args_o, 'UniformOutput', false);

        obj.func = casadi.Function(name, sym_i, sym_o, n_i, n_o, varargin{:});
    end

    function nm = get.name(obj)
        % Return function name.
        nm = obj.func.name;
    end

    function argout = call(obj, argin, argout)
        % Evaluate casadi function object.
        in = structfun(@get_param, argin, 'UniformOutput', false);

        % call casadi function
        out = call(obj.func, in);

        % return result
        for fn = fieldnames(out)
            argout.(fn{:}) = set_param(argout.(fn{:}), out.(fn{:}));
        end
    end
end

end
