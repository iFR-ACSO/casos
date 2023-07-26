classdef (Abstract) FunctionConstArgument < casos.package.functions.FunctionArgument
% Function argument with constant type (DM, SX, MX).

properties (SetAccess=protected)
    value = [];
end

methods
    function arg = FunctionConstArgument(varargin)
        % New argument.
        arg@casos.package.functions.FunctionArgument(varargin{:});
    end

    %% Interface
    % obj = set_value(obj,val);

    function tf = has_value(obj)
        % Check if argument has value.
        tf = ~isempty(obj.value);
    end

    function tf = is_symbolic(obj)
        % Check if argument is symbolic variable.
        tf = is_symbolic(obj.expr);
    end

    function val = get_value(obj)
        % Return value of argument.
        val = obj.value;
    end

    function par = get_param(obj)
        % Return actual parameter.
        par = get_value(obj);
    end

    function obj = set_param(obj,par)
        % Set actual parameter.
        obj.value = casadi.DM(par); % always return constant argument as DM
    end

    function sym = get_symbolic(obj)
        % Return formal parameter.
        sym = obj.expr;
    end
end

end
