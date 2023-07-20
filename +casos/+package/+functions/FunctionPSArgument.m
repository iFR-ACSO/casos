classdef FunctionPSArgument < casos.package.functions.FunctionArgument
% Function argument with polynomial (PS) type.

properties (Constant)
    type = casos.package.functions.FunctionArgumentType.PS;
end

properties (SetAccess=protected)
    value = [];
end

properties (Access=private)
    coeffs;
    monoms;
end

methods
    function arg = FunctionPSArgument(varargin)
        % New argument.
        arg@casos.package.functions.FunctionArgument(varargin{:});

        % get coefficients and monomials of expression
        [arg.coeffs,arg.monoms] = poly2basis(arg.expr);
    end

    %% Interface
    function tf = has_value(obj)
        % Check if argument has value.
        tf = ~isempty(obj.value);
    end

    function tf = is_symbolic(obj)
        % Check if argument is symbolic variable.
        tf = is_symbolic(obj.coeffs);
    end

    function val = get_value(obj)
        % Return value of argument.
        val = obj.value;
    end

    function obj = set_value(obj,val)
        % Set value of argument.
        obj.value = casos.PS(val);
    end

    function par = get_param(obj)
        % Return actual parameter.
        par = poly2basis(obj.value,obj.monoms);
    end

    function obj = set_param(obj,par)
        % Set actual parameter.
        obj.value = reshape(par,1,length(obj.monoms))*obj.monoms;
    end

    function sym = get_symbolic(obj)
        % Return formal parameter.
        sym = obj.coeffs;
    end
end

end
