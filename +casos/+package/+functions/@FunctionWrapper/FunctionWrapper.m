classdef (Abstract) FunctionWrapper
% Wrap a function interface.
    
properties (Dependent)
    class_name;

    name;

    % input arguments
    n_in;
%     name_in;
%     monomials_in;
%     default_in;
%     size_in;

    % output arguments
    n_out;
%     name_out;
%     monomials_out;
%     size_out;
    stats;
end

properties (Access = private)
    wrap;
end

methods (Access=protected,Static)
    type = parse_argument(expr);
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

    function n = get.n_in(obj)
        % Return number of inputs.
        n = get_n_in(obj.wrap);
    end

    function nm = name_in(obj,varargin)
        % Return name of inputs.
        nm = get_name_in(obj.wrap,varargin{:});
    end

    function z = monomials_in(obj,i)
        % Return monomials of inputs.
        z = get_monomials_in(obj.wrap,i);
    end

    function val = default_in(obj,i)
        % Return default values for inputs.
        val = get_default_in(obj.wrap,i);
    end

    function sz = size_in(obj,i)
        % Return size of inputs.
        sz = get_size_in(obj.wrap,i);
    end

    function i = index_in(obj,str)
        % Return index of inputs.
        i = get_index_in(obj.wrap,str);
    end

    function n = get.n_out(obj)
        % Return number of outputs.
        n = get_n_out(obj.wrap);
    end

    function nm = name_out(obj,varargin)
        % Return name of outputs.
        nm = get_name_out(obj.wrap,varargin{:});
    end

    function z = monomials_out(obj,i)
        % Return monomials of outputs.
        z = get_monomials_out(obj.wrap,i);
    end

    function sz = size_out(obj,i)
        % Return size of outputs.
        sz = get_size_out(obj.wrap,i);
    end

    function i = index_out(obj,str)
        % Return index of outputs.
        i = get_index_out(obj.wrap,str);
    end

    function s = get.stats(obj)
        % Return stats.
        s = get_stats(obj.wrap);
    end

    function print_options(obj)
        % Print list of options.
        print_options(obj.wrap);
    end

    function print_option(obj,name)
        % Print information about an option.
        print_option(obj.wrap,name);
    end

    function tf = has_option(obj,name)
        % Check if option "name" exists.
        tf = has_option(obj.wrap,name);
    end

    function out = call(obj,args)
        % Evaluate function for given arguments.
        if iscell(args)
            assert(length(args) == obj.n_in, 'Incorrect number of inputs: Expected %d, got %d.', obj.n_in, length(args));

            out = call(obj.wrap,args);
            return
        end

        % else
        assert(isstruct(args),'Arguments must be given as struct.')

        argin = cell(obj.n_in,1);

        % name of arguments
        fn_arg = fieldnames(args);

        % index of arguments
        idx_in = cellfun(@(arg) obj.index_in(arg), fn_arg);
        
        % find arguments
        L = ismember(0:obj.n_in-1, idx_in);

        % default values
        argin(~L) = arrayfun(@(i) obj.default_in(i), find(~L)-1, 'UniformOutput', false);

        % assign arguments
        argin(idx_in+1) = struct2cell(args);
        
        % call function with cell
        argout = call(obj.wrap,argin);

        % name of outputs
        fn_out = arrayfun(@(i) obj.name_out(i), 0:obj.n_out-1, 'UniformOutput', false);

        % parse outputs
        out = cell2struct(argout(:),fn_out(:));
    end

    function varargout = subsref(obj,L)
        % Subscripted reference.
        if length(L) > 1 || ~strcmp(L.type,'()')
            % fall back to builtin function
            [varargout{1:nargout}] = builtin('subsref',obj,L);

            return
        end

        % else:
        if isempty(L.subs)
            % empty function call
            args = struct;
            
        elseif ischar(L.subs{1})
            % name-value pairs
            fn_i = L.subs(1:2:end);
            L.subs(1:2:end) = [];

            assert(length(fn_i) == length(L.subs), 'Name-value syntax requires same number of names and values.')

            % assign inputs
            args = cell2struct(L.subs(:),fn_i(1:length(L.subs)));
        else
            % multiple inputs
            args = L.subs;
        end

        % call function
        out = call(obj,args);

        if isstruct(out)
            varargout = {out};
        else
            % return multiple outputs
            varargout = out;
        end
    end
end

methods (Access={?casos.package.functions.FunctionInterface})
    %% Friend interface
    function f = substitute(obj,varargin)
        % Substitute variables.
        f = obj;
        f.wrap = substitute(obj.wrap,varargin{:});
    end
end

end
