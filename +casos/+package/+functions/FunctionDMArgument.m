classdef FunctionDMArgument < casos.package.functions.FunctionConstArgument
% Function argument of type DM (or double).

properties (Constant)
    type = casos.package.functions.FunctionArgumentType.DM;
end

methods
    function arg = FunctionDMArgument(varargin)
        % New argument.
        arg@casos.package.functions.FunctionConstArgument(varargin{:});
    end

    %% Interface
    function obj = set_value(obj,val)
        % Set value of argument.
        obj.value = casadi.DM(val);
    end
end

end
