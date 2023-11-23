classdef Options
% Options.

properties (Access=private)
    entries = table({},'VariableNames',{'description'});
end

methods
    function obj = Options(args)
        % Create new Options instance from cell or struct.
        if nargin < 1
            % nothing to do
        elseif isa(args,'casos.package.Options')
            % copy constructor
            obj = args;
        elseif iscell(args) && size(args,2) == 2
            % cell array
            obj.entries = table(args(:,2),'VariableNames',{'description'},'RowNames',args(:,1));
        else
            error('Options undefined for input of type %s.',class(args))
        end
    end

    function obj = vertcat(obj,args)
        % Concatenate Options.
        args = casos.package.Options(args);

        % concatenate fields
        obj.entries = [obj.entries; args.entries];
    end

    function tf = has(obj,name)
        % Check if option "name" exists.
        tf = ismember(name,obj.entries.Row);
    end

    function str = info(obj,name)
        % Return description for option
        assert(has(obj,name), 'Unkown option "%s".',name);
        % description
        str = obj.entries(name,:).description{:};
    end

    function check(obj,opts)
        % Check option structure.
        fn_opts = fieldnames(opts);
        % check if options exist
        tf = cellfun(@(fn) has(obj,fn), fn_opts);
        % throw error
        assert(all(tf), 'Unkown option "%s".', fn_opts{find(tf,1)});
    end

    function print_all(obj)
        % Print all options.
        disp(obj.entries)
    end

    function print_one(obj,name)
        % Print one option.
        assert(has(obj,name), 'Unkown option "%s".',name)
        % disply single option
        disp(obj.entries(name,:))
    end

    function disp(obj)
        % Display Options.
        disp(obj.entries)
    end
end

end
