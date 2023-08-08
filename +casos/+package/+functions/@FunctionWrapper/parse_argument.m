function [type,args] = parse_argument(name,expr)
% Return type and argument object for expression.

if isdouble(expr) || isa(expr,'casadi.DM')
    args = casos.package.functions.FunctionDMArgument(name,expr);

elseif isa(expr,'casadi.SX')
    args = casos.package.functions.FunctionSXArgument(name,expr);

elseif isa(expr,'casadi.MX')
    args = casos.package.functions.FunctionMXArgument(name,expr);

elseif isa(expr,'casos.PS')
    args = casos.package.functions.FunctionPSArgument(name,expr);

else
    error('Function undefined for class %s of input (%s).',class(expr),name);
end

type = args.type;

end