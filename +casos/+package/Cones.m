classdef Cones
% Cones.

properties (Access=private)
    entries = table({},{},'VariableNames',{'type' 'description'});
end

properties (Constant)
    % Predefined cone types
    LIN = {'lin' 'NUM' 'Linear inequalities (element-wise).'};
    LOR = {'lor' 'LIST' 'Lorentz (quadratic, second-order) cone.'};
    ROT = {'rot' 'LIST' 'Rotated Lorentz cone.'};
    PSD = {'psd' 'LIST' 'Cone of positive semidefinite matrices.'};
    % Polynomial cones
    SOS = {'sos' 'NUM' 'Cone of sum-of-squares polynomials.'}
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

        % legacy
        switch (name)
            case 'l', warning('Cone "l" is not supported, use "lin".')
            case 'q', warning('Cone "q" is not supported, use "lor".')
            case 'r', warning('Cone "r" is not supported, use "rot".')
            case 's', warning('Cone "s" is not supported, use "psd" or "sos".')
        end
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

    function l = get_length(obj,K,name,no_check)
        % Return length of cones.
        if nargin < 3
            % return total length
            l = sum(cellfun(@(fn) get_length(obj,K,fn,true), fieldnames(K)));
            return
        end

        % else:
        % return length of specified cone
        assert(has(obj,name) || (nargin > 3 && no_check), 'Unknown cone "%s".',name)
        
        % check if cone is specified
        if ~isfield(K,name)
            l = 0;
            return
        end

        % else
        switch (name)
            case 'psd', l = sum(K.(name).^2);
            otherwise,  l = sum(K.(name));
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
