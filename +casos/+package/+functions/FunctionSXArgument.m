classdef FunctionSXArgument < casos.package.functions.FunctionConstArgument
% Function argument of type SX.

properties (Constant)
    type = casos.package.functions.FunctionArgumentType.SX;
end

methods
    function arg = FunctionSXArgument(varargin)
        % New argument.
        arg@casos.package.functions.FunctionConstArgument(varargin{:});
    end

    %% Interface
    function obj = set_value(obj,val)
        % Set value of argument.
        obj.value = casadi.SX(val);
    end
end

end
