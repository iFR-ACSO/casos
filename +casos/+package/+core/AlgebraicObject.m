classdef (Abstract, InferiorClasses = {?casadi.DM, ?casadi.SX, ?casadi.MX}) ...
        AlgebraicObject < matlab.mixin.indexing.RedefinesParen
% Base class for algebraic objects.

methods (Abstract)
    % check if object represents indeterminate variables
    tf = is_indet(obj);
end

methods
    %% Algebraic interface
    function t = cat(N,varargin)
        % Concatenate algebraic objects.
        args = cellfun(@(v) casos.package.polynomial(v), varargin, 'UniformOutput',false);
        t = cat(N,args{:});
    end
end

methods (Access=protected)
    %% protected RedefinesParen interface
    function n = parenListLength(obj,indexOp,context)
        % Return length of indexing operation into algebraic object.
        if length(indexOp) < 2, n = 1; else, n = listLength(obj,indexOp(2:end),context); end
    end
end

end