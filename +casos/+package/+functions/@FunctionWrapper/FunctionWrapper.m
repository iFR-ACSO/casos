classdef (Abstract) FunctionWrapper
% Wrap a function interface.
    
properties (Dependent)
    class_name;
    name;
end

properties (Access = private)
    wrap;
end

methods (Access=protected,Static)
    [type,args] = parse_argument(name,expr);
end

methods
    function f = FunctionWrapper(wrap)
        % Superclass constructor.
        f.wrap = wrap;
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
        argin = set_value_in(obj.wrap,args);

        % call wrapped function
        argout = call(obj.wrap,argin);

        % parse outputs
        out = get_value_out(obj.wrap,argout);
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

            fn_i = get_name_in(obj.wrap);
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
