classdef (Abstract) FunctionInterface < casos.package.functions.FunctionCommon
% Interface for callable functions.

properties (Abstract,SetAccess=protected)
    class_name;
end

properties (SetAccess = private)
    name;
end

methods (Abstract)
    argout = call(obj,argin);

    % inputs
    n = get_n_in(obj);
    s = get_name_in(obj,i);
    z = get_sparsity_in(obj,i);
%    z = get_monomials_in(obj,i);
%    v = get_default_in(obj,i);
%    s = get_size_in(obj,i);
    i = get_index_in(obj,str);
    % outputs
    n = get_n_out(obj);
    z = get_sparsity_out(obj,i);
    s = get_name_out(obj,i);
%    z = get_monomials_out(obj,i);
%    s = get_size_out(obj,i);
    i = get_index_out(obj,str);
end

methods
    function obj = FunctionInterface(name,varargin)
        % Superclass constructor.
        obj@casos.package.functions.FunctionCommon(varargin{:});
        obj.name = name;
    end

    function z = get_monomials_in(obj,i)
        % Return monomials of input arguments.
        warning('DEPRECATED. Use sparsity instead.')
        z = monomials(get_sparsity_in(obj,i));
    end

    function val = get_default_in(obj,i) %#ok<INUSD>
        % Return default values of input arguments.
        val = 0;
    end

    function sz = get_size_in(obj,i)
        % Return size of input arguments.
        sz = size(get_sparsity_in(obj,i));
    end

    function z = get_monomials_out(obj,i)
        % Return monomials of output arguments.
        warning('DEPRECATED. Use sparsity instead.')
        z = monomials(get_sparsity_out(obj,i));
    end

    function sz = get_size_out(obj,i)
        % Return size of output arguments.
        sz = size(get_sparsity_out(obj,i));
    end

    function s = get_stats(~)
        % Return empty stats.
        s = struct;
    end
end

methods (Access={?casos.package.functions.FunctionCommon, ?casos.package.functions.FunctionWrapper})
    %% Friend interface
    function substitute(obj,varargin)
        % Substitute variables.
        error('Method SUBSTITUTE not supported for %s.',obj.class_name)
    end
end

end
