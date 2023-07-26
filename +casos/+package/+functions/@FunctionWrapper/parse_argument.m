function [type,args] = parse_argument(expr,name)
% Return type and argument object for expression.

if isdouble(expr) || isa(expr,'casadi.DM')
    args = casos.package.functions.FunctionDMArgument(expr,name);

elseif isa(expr,'casadi.SX')
    args = casos.package.functions.FunctionSXArgument(expr,name);

elseif isa(expr,'casadi.MX')
    args = casos.package.functions.FunctionMXArgument(expr,name);

elseif isa(expr,'casos.PS')
    args = casos.package.functions.FunctionPSArgument(expr,name);

else
    error('Function undefined for class %s of input (%s).',class(expr),name);
end

type = args.type;

end