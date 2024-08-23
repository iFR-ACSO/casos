classdef FunctionArgumentType
% Possible types of arguments for functions.

    properties
        classname;
    end

    enumeration
        % vectors & matrices
        DM ('casadi.DM')
        SX ('casadi.SX')
        MX ('casadi.MX') 

        % polynomials
        PD ('casos.PD')
        PS ('casos.PS')

        % operators
        PDOperator ('casos.PDOperator')
        PSOperator ('casos.PSOperator')
    end

    methods
        function type = FunctionArgumentType(cls)
            % Create new function argument type.
            type.classname = cls;
        end

        function tf = isOfType(type,var)
            % Check if variable is of given type.
            tf = isa(var, type.classname);
        end
    end
end
