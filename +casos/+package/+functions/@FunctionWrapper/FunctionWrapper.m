classdef (Abstract) FunctionWrapper
% Wrap a function interface.
    
properties (Dependent)
    class_name;
    name;
end

properties (SetAccess = private)
    arg_i;
    arg_o;
end

properties (Access = private)
    wrap;
end

methods (Access=protected,Static)
    [type,args] = parse_argument(expr,name);
end

methods
    function f = FunctionWrapper(wrap,arg_i,arg_o,name_i,name_o)
        % Superclass constructor.
        f.wrap = wrap;
        
        % input/output arguments
        f.arg_i = cell2struct(arg_i,name_i);
        f.arg_o = cell2struct(arg_o,name_o);
    end

    function cls = get.class_name(obj)
        % Return function class.
        cls = obj.wrap.class_name;
    end

    function nm = get.name(obj)
        % Return function name.
        nm = obj.wrap.name;
    end

    function out = call(obj,args)
        % Evaluate function for given arguments.
        assert(isstruct(args),'Arguments must be given as struct.')
        
        % parse inputs by name
        argin = obj.arg_i;
        for fn = fieldnames(args)
            assert(ismember(fn,fieldnames(argin)), 'Unexpected input (%s)', fn{:})

            argin.(fn{:}) = set_value(argin.(fn{:}), args.(fn{:}));
        end

        % call wrapped function
        argout = call(obj.wrap,argin,obj.arg_o);

        % parse outputs
        out = [];
        for fn = fieldnames(argout)
            assert(has_value(argout.(fn{:})), 'Output (%s) not assigned after function call.', fn{:})

            out.(fn{:}) = get_value(argout.(fn{:}));
        end
    end

    function varargout = subsref(obj,L)
        % Subscripted reference.
        if length(L) > 1 || ~strcmp(L.type,'()')
            % fall back to builtin function
            [varargout{1:nargout}] = builtin('subsref',obj,L);

            return
        end

        % else:
        fn_i = fieldnames(obj.arg_i);
        fn_o = fieldnames(obj.arg_o);

        % assign inputs
        args = cell2struct(L.subs,fn_i(1:length(L.subs)));

        % call function
        out = call(obj,args);

        % return first output
        varargout = {out.(fn_o{1})};
    end
end

end
