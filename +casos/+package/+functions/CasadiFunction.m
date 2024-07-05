classdef CasadiFunction < casos.package.functions.FunctionInternal
% Casadi function interface.

properties (Access=private)
    func;
end

properties (Dependent,SetAccess=private)
    class_name;
end

properties (Constant,Access=protected)
    allow_eval_on_basis = false;
end

methods
    function obj = CasadiFunction(name, ex_i, ex_o, name_i, name_o, varargin)
        % Create new casadi function object.
        if isa(name,'casadi.Function')
            % existing casadi Function object
            func = name;
            name = func.name;
        else
            % create new casadi Function
            func = casadi.Function(name, ex_i, ex_o, name_i, name_o, varargin{:});
        end

        obj@casos.package.functions.FunctionInternal(name);

        obj.func = func;
    end

    function cls = get.class_name(obj)
        % Return function class name.
        cls = obj.func.class_name;
    end

    %% Options & Cones
    function print_options(obj)
        % Print list of options.
        print_options(obj.func);
    end

    function print_option(obj,name)
        % Print information about an option.
        print_option(obj.func,name);
    end

    function tf = has_option(obj,name)
        % Check if option "name" exists.
        tf = has_option(obj.func,name);
    end

    function print_cones(obj)
        % Print list of supported cones.
        if isa(obj.func,'casos.package.functions.FunctionCommon')
            print_cones(obj.func);
        else
            error('Not implemented for class %s.',obj.class_name)
        end
    end

    function print_cone(obj,name)
        % Print information about a cone.
        if isa(obj.func,'casos.package.functions.FunctionCommon')
            print_cone(obj.func,name);
        else
            error('Not implemented for class %s.',obj.class_name)
        end
    end

    function tf = has_cone(obj,name)
        % Check if cone "name" is supported.
        if isa(obj.func,'casos.package.functions.FunctionCommon')
            tf = has_cone(obj.func,name);
        else
            error('Not implemented for class %s.',obj.class_name)
        end
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

    function z = get_sparsity_in(obj,i)
        % All inputs are of degree zero.
        z = casos.Sparsity(sparsity_in(obj.func,i));
    end

    function val = get_default_in(obj,i)
        % Default inputs.
        val = default_in(obj.func,i);
    end

    function i = get_index_in(obj,str)
        % Index of inputs.
        i = index_in(obj.func,str);
    end

    function n = get_n_out(obj)
        % Number of outputs.
        n = n_out(obj.func);
    end

    function str = get_name_out(obj,varargin)
        % Name of outputs.
        str = name_out(obj.func,varargin{:});
    end

    function z = get_sparsity_out(obj,i)
        % All outputs are of degree zero.
        z = casos.Sparsity(sparsity_out(obj.func,i));
    end

    function i = get_index_out(obj,str)
        % Index of outputs.
        i = index_out(obj.func,str);
    end

    function argout = call(obj, argin)
        % Evaluate casadi function object.
        argout = call(obj.func, argin);
    end

    function s = get_stats(obj)
        % Return stats.
        s = stats(obj.func);
    end
end

end
