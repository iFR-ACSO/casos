classdef FunctionMXArgument < casos.package.functions.FunctionConstArgument
% Function argument of type MX.

properties (Constant)
    type = casos.package.functions.FunctionArgumentType.MX;
end

methods
    function arg = FunctionMXArgument(varargin)
        % New argument.
        arg@casos.package.functions.FunctionConstArgument(varargin{:});
    end

    %% Interface
    function obj = set_value(obj,val)
        % Set value of argument.
        obj.value = casadi.MX(val);
    end
end

end
