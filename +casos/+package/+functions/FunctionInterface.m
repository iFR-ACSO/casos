classdef (Abstract) FunctionInterface
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
    z = get_monomials_in(obj,i);
%     v = get_default_in(obj,i);
    s = get_size_in(obj,i);
    % outputs
    n = get_n_out(obj);
    s = get_name_out(obj,i);
    z = get_monomials_out(obj,i);
    s = get_size_out(obj,i);
end

methods
    function obj = FunctionInterface(name)
        % Superclass constructor.
        obj.name = name;
    end

    function val = get_default_in(obj,i)
        % Return default values of input arguments.
        val = 0;
    end
end

end
