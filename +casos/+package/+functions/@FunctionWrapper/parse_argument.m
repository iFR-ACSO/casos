function type = parse_argument(expr)
% Return type of expression.

if isnumeric(expr) || isa(expr,'casadi.DM')
    type = casos.package.functions.FunctionArgumentType.DM;

elseif isa(expr,'casadi.SX')
    type = casos.package.functions.FunctionArgumentType.SX;

elseif isa(expr,'casadi.MX')
    type = casos.package.functions.FunctionArgumentType.MX;

elseif isa(expr,'casos.PD')
    type = casos.package.functions.FunctionArgumentType.PD;

elseif isa(expr,'casos.PS')
    type = casos.package.functions.FunctionArgumentType.PS;

elseif isa(expr,'casos.PDOperator')
    type = casos.package.functions.FunctionArgumentType.PDOperator;

elseif isa(expr,'casos.PSOperator')
    type = casos.package.functions.FunctionArgumentType.PSOperator;

else
    error('Function undefined for class %s of input (%s).',class(expr),name);
end

end