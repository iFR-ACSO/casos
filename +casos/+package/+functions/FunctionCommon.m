classdef (Abstract) FunctionCommon < handle
% Common superclass for casos and casadi-derived functions.

properties (Constant, Access=protected)
    options = casos.package.Options({'error_on_fail', 'Throw exceptions when function evaluation fails (default true).'});
end

properties (Access=protected)
    opts = struct;
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.functions.FunctionCommon.options;
    end
end

methods
    function obj = FunctionCommon(opts)
        % Instantiate base class.
        if nargin > 0
            % assert that all options exist
            check(obj.get_options, opts);
            % store options
            obj.opts = opts;
        end

        % default options
        if ~isfield(obj.opts,'error_on_fail')
            obj.opts.error_on_fail = true;
        end
    end

    function print_options(obj)
        % Print list of options.
        print_all(obj.get_options);
    end

    function print_option(obj,name)
        % Print information about an option.
        print_one(obj.get_options,name);
    end

    function tf = has_option(obj,name)
        % Check if option "name" exists.
        tf = has(obj.get_options,name);
    end
end

methods (Abstract, Access={?casos.package.functions.FunctionCommon, ?casos.package.functions.FunctionWrapper})
    %% Friend interface
    f = substitute(obj,varargin);
end

end
