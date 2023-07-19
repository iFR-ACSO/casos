classdef (Abstract) FunctionArgument
% Base class of function arguments.

properties (Abstract,Constant)
    type;
end

properties (SetAccess=private)
    expr;
    name;
end

properties (Dependent)
    signature;
end

methods (Abstract)
    tf = has_value(obj);
    tf = is_symbolic(obj);

    val = get_value(obj);
    obj = set_value(obj,val);

    par = get_param(obj);
    obj = set_param(obj,par);

    sym = get_symbolic(obj);
end

methods
    function arg = FunctionArgument(expr,name)
        % Superclass constructor.
        arg.expr = expr;
        arg.name = name;
    end

    function sg = get.signature(obj)
        % Return signature.
        sg = sprintf('%s:%s', obj.name, obj.type);
    end
end

end
