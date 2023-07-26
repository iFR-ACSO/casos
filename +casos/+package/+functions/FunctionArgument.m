classdef (Abstract) FunctionArgument
% Base class of function arguments.

properties (Abstract,Constant)
    type;
end

properties (SetAccess=private)
    expr;
    name;
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

    function sg = get_signature(obj)
        % Return signature.
        dim = get_dimensions(obj);

        if ~isempty(dim)
            dim = {sprintf('[%s]',[dim{:}])};
        end

        sg = sprintf('%s', obj.name, [dim{:}]);
    end
end

methods (Access=protected)
    function dim = get_dimensions(obj)
        % Return dimensions of argument.
        [n,m] = size(obj.expr);

        if n == 1 && m == 1
            % don't show dimensions
            dim = {};
        elseif m == 1
            % show number of rows
            dim = compose('%d',n);
        elseif n > 0 || m > 0
            % show size
            dim = compose('%dx%d',n,m);
        else
            % show empty dimensions
            dim = {''};
        end
    end
end

end
