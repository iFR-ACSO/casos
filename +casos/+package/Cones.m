classdef Cones
% Cones.

properties (Access=private)
    entries = table({},{},'VariableNames',{'type' 'description'});
end

methods
    function obj = Cones(args)
        % Create new Cones instance from cell.
        if nargin < 1
            % nothing to do
        elseif isa(args,'casos.package.Cones')
            % copy constructor
            obj = args;
        elseif iscell(args) && size(args,2) == 3
            % cell array
            obj.entries = table(args(:,2),args(:,3),'VariableNames',{'type' 'description'},'RowNames',args(:,1));
        else
            error('Cones undefined for input of type %s.',class(args))
        end
    end

    function obj = vertcat(obj,args)
        % Concatenate Cones.
        args = casos.package.Cones(args);

        % concatenate fields
        obj.entries = [obj.entries; args.entries];
    end

    function tf = has(obj,name)
        % Check if cone "name" exists.
        tf = ismember(name,obj.entries.Row);
    end

    function str = info(obj,name)
        % Return description for cone
        assert(has(obj,name), 'Unknown cone "%s".',name);
        % description
        str = obj.entries(name,:).description{:};
    end

    function check(obj,opts)
        % Check cone structure.
        fn_cone = fieldnames(opts);
        % check if cones exist
        tf = cellfun(@(fn) has(obj,fn), fn_cone);
        % throw error
        assert(all(tf), 'Unknown cone "%s".', fn_cone{find(~tf,1)});
    end

    function l = get_length(obj,K,name)
        % Return length of specified cone.
        assert(has(obj,name), 'Unknown cone "%s".',name)
        % check if cone is specified
        if isfield(K,name)
            l = sum(K.(name));
        else
            l = 0;
        end
    end

    function d = get_dimension(obj,K,name)
        % Return dimension of specified cone.
        assert(has(obj,name), 'Unknown cone "%s".',name)
        % check if cone is specified
        if isfield(K,name)
            d = K.(name);
        elseif strcmpi(obj.entries(name,:).type,'list')
            d = [];
        else
            d = 0;
        end
    end

    function print_all(obj)
        % Print all cones.
        disp(obj.entries)
    end

    function print_one(obj,name)
        % Print one cone.
        assert(has(obj,name), 'Unknown cone "%s".',name)
        % disply single cone
        disp(obj.entries(name,:))
    end

    function disp(obj)
        % Display cones.
        disp(obj.entries)
    end
end

end
