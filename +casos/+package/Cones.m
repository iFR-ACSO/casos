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
    DD  = {'dd'  'LIST' 'Cone of symmetric diagonally dominat matrices'};
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
        str = obj.entries.description{name};
    end

    function check(obj,K)
        % Check cone structure.
        fn_cone = fieldnames(K);
        % check if cones exist
        cellfun(@(fn) check_cone(obj,K,fn), fn_cone);
    end

    function check_cone(obj,K,name)
        % Check specified cone.
        assert(has(obj,name), 'Unknown cone "%s".', name);

        if ~isfield(K,name)
            % nothing to do
            return
        end

        % else
        switch lower(obj.entries.type{name})
            case 'num'
                assert(isscalar(K.(name)), 'Cone "%s" must be specified as scalar.', name);
            case 'list'
                assert(isvector(K.(name)), 'Cone "%s" must be specified as vector.', name);
            case {'struct' 'cones'}
                assert(isstruct(K.(name)), 'Cone "%s" must be specified as struct.', name);
        end
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
            case {'psd' 'dd'}, l = sum(K.(name).^2);
            otherwise,  l = sum(K.(name));
        end
    end

    function d = get_dimension(obj,K,name)
        % Return dimension of specified cone.
        assert(has(obj,name), 'Unknown cone "%s".',name)
        % check if cone is specified
        if isfield(K,name)
            d = K.(name);
        elseif strcmpi(obj.entries.type{name},'list')
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
