classdef (Abstract) FunctionInterface
% Interface for callable functions.

properties (Abstract,SetAccess=protected)
    class_name;
end

methods (Abstract)
    argout = call(obj,argin,argout);
end

end
