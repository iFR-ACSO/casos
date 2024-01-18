classdef (Abstract) SolverCommon < casos.package.functions.FunctionCommon
% Common solver interface.

properties (Constant,Access=protected)
    solver_options = [casos.package.functions.FunctionInterface.options
        {'Kx', 'Cone description for state constraints.'
         'Kc', 'Cone description for constraint function.'}
    ];
end

methods (Static,Abstract)
    cones = get_cones;
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SolverCommon.solver_options;
    end
end

methods
    function obj = SolverCommon(varargin)
        % Instantiate base class.
        obj@casos.package.functions.FunctionCommon(varargin{:});

        if isfield(obj.opts,'Kx')
            check(obj.get_cones, obj.opts.Kx);
        end
        if isfield(obj.opts,'Kc')
            check(obj.get_cones, obj.opts.Kc);
        end
    end

    function print_cones(obj)
        % Print list of supported cones.
        print_all(obj.get_cones);
    end

    function print_cone(obj,name)
        % Print information about a cone.
        print_one(obj.get_cones,name);
    end

    function tf = has_cone(obj,name)
        % Check if cone "name" is supported.
        tf = has(obj.get_cones,name);
    end
end

end
