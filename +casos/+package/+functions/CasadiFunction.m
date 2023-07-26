classdef CasadiFunction < casos.package.functions.FunctionInterface
% Casadi function interface.

properties (Access=private)
    func;
end

properties (Dependent)
    class_name;
    name;
end

methods
    function obj = CasadiFunction(name, args_i, args_o, varargin)
        % Create new casadi function object.
        [ex_i,n_i] = cellfun(@(arg) deal(arg.expr,arg.name), args_i, 'UniformOutput', false);
        [ex_o,n_o] = cellfun(@(arg) deal(arg.expr,arg.name), args_o, 'UniformOutput', false);

        obj.func = casadi.Function(name, ex_i, ex_o, n_i, n_o, varargin{:});
    end

    function cls = get.class_name(obj)
        % Return function class name.
        cls = obj.func.class_name;
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
