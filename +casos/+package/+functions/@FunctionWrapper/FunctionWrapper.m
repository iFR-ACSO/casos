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
    [type,args] = parse_argument(name,expr);
end

methods
    function f = FunctionWrapper(wrap,arg_i,arg_o,name_i,name_o)
        % Superclass constructor.
        f.wrap = wrap;
        
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
        if ischar(L.subs{1})
            % name-value pairs
            namevalue = true;

            fn_i = L.subs(1:2:end);
            L.subs(1:2:end) = [];

            assert(length(fn_i) == length(L.subs), 'Name-value syntax requires same number of names and values.')
        else
            namevalue = false;

            fn_i = fieldnames(obj.arg_i);
        end
        
        % assign inputs
        args = cell2struct(L.subs,fn_i(1:length(L.subs)));

        % call function
        out = call(obj,args);

        if namevalue
            varargout = {out};
        else
            % return first output
            varargout = struct2cell(out);
        end
    end
end

end
