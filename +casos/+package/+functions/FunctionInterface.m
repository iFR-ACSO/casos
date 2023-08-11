classdef (Abstract) FunctionInterface
% Interface for callable functions.

properties (Abstract,SetAccess=protected)
    class_name;
    name;
end

properties (SetAccess = private)
    arg_i;
    arg_o;
end

methods (Abstract)
    argout = call(obj,argin);
end

methods
    function f = FunctionInterface(arg_i,arg_o,name_i,name_o)
        % Superclass constructor.

        % input/output arguments
        if nargin < 4
            assert(istruct(arg_i) && istruct(arg_o), 'Arguments must be structures.')
            f.arg_i = arg_i;
            f.arg_o = arg_o;
        else
            f.arg_i = cell2struct(arg_i,name_i);
            f.arg_o = cell2struct(arg_o,name_o);
        end
    end

    function str = get_name_in(obj)
        % Return name of input arguments.
        str = fieldnames(obj.arg_i);
    end

    function str = get_name_out(obj)
        % Return name of output arguments.
        str = fieldnames(obj.arg_o);
    end

    function arg = set_value_in(obj,values)
        % Set input values before call.
        arg = obj.arg_i;
        for fn = fieldnames(values)
            assert(ismember(fn,fieldnames(arg)), 'Unexpected input (%s)', fn{:})

            arg.(fn{:}) = set_value(arg.(fn{:}), values.(fn{:}));
        end
    end

    function values = get_value_out(~,arg)
        % Get output values after call.
        values = [];
        for fn = fieldnames(arg)
            assert(has_value(arg.(fn{:})), 'Output (%s) not assigned after function call.', fn{:})

            values.(fn{:}) = get_value(arg.(fn{:}));
        end
    end
end

end