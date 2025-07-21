classdef (Abstract) SolverCallback < casadi.Callback & casos.package.solvers.SolverCommon
% Common Callback interface for custom solvers.

properties (Abstract, Access=protected)
    fhan;
    ghan;
end

properties (Constant, Access=protected)
%     solver_options = casos.package.functions.FunctionCommon.options;
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SolverCallback.solver_options;
    end
end

methods
    function obj = SolverCallback(varargin)
        obj@casos.package.solvers.SolverCommon(varargin{:});
    end

    %% Common function interface
    function print_options(obj)
        % Print list of options.
        print_options@casos.package.functions.FunctionCommon(obj);
    end

    function print_option(obj,name)
        % Print information about an option.
        print_option@casos.package.functions.FunctionCommon(obj,name);
    end

    function tf = has_option(obj,name)
        % Check if option "name" exists.
        tf = has_option@casos.package.functions.FunctionCommon(obj,name);
    end

    %% Common Callback interface
    function n = get_n_in(obj)
        % Return number of inputs.
        n = n_in(obj.fhan);
    end

    function n = get_n_out(obj)
        % Return number of outputs.
        n = n_out(obj.ghan);
    end

    function str = get_name_in(obj,i)
        % Return names of input arguments.
        str = name_in(obj.fhan,i);
    end

    function str = get_name_out(obj,i)
        % Return names of output arguments.
        str = name_out(obj.ghan,i);
    end

    function sp = get_sparsity_in(obj,i)
        % Return sparsity of input arguments.
        sp = sparsity_in(obj.fhan,i);
    end

    function sp = get_sparsity_out(obj,i)
        % Return sparsity of output arguments.
        sp = sparsity_out(obj.ghan,i);
    end
end

end
