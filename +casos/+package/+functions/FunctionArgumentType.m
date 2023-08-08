classdef FunctionArgumentType
% Possible types of arguments for functions.

    properties
        classname;
        construct;
    end

    enumeration
        DM ('casadi.DM', @casos.package.functions.FunctionDMArgument)
        SX ('casadi.SX', @casos.package.functions.FunctionSXArgument)
        MX ('casadi.MX', @casos.package.functions.FunctionMXArgument) 
        PS ('casos.PS', @casos.package.functions.FunctionPSArgument)
    end

    methods
        function type = FunctionArgumentType(cls,cst)
            % Create new function argument type.
            type.classname = cls;
            type.construct = cst;
        end

        function tf = isOfType(type,var)
            % Check if variable is of given type.
            tf = isa(var, type.classname);
        end

        function arg = newArgument(type,varargin)
            % Create new function argument of given type.
            arg = feval(type.construct, varargin{:});
        end
    end
end
