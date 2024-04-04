classdef (Abstract) SolverCommon < casos.package.functions.FunctionCommon
% Common solver interface.

properties (Constant,Access=protected)
    solver_options = [casos.package.functions.FunctionCommon.options
        {'Kx', 'Cone description for state constraints.'
         'Kc', 'Cone description for constraint function.'
         'add_lbc', 'something 1'
         'add_ubc', 'something 2'}
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
        disp('Supported Cones:')
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

    function print_options(obj)
        % Overwrite FunctionCommon#print_options
        print_options@casos.package.functions.FunctionCommon(obj);

        print_cones(obj);
    end

    function print_option(obj,name)
        % Overwrite FunctionCommon#print_option
        names = split(name,'.');
        % print specified option
        print_option@casos.package.functions.FunctionCommon(obj,names{1});

        % print specified cone, if any
        if ~ismember(names{1},{'Kx' 'Kc'}) %TODO: check for type rather than name
            assert(length(names) == 1,'Undefined behaviour.')
        elseif length(names) > 1
            print_cone(obj,names{2});
        else
            print_cones(obj);
        end
    end
end

end
