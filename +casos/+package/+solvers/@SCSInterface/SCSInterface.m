classdef (Sealed) SCSInterface < casos.package.solvers.ConicSolver
% Interface for conic solver SCS.

properties (Access=protected)
    fhan;
    ghan;
    cone;
end

properties (Access=private)
    info;
end

methods
%     init(obj);
    argout = eval(obj,argin);

    function obj = SCSInterface(name,conic,varargin)
        % Construct SCS interface.
        obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});
    end

    function s = stats(obj)
        % Return stats.
        s = obj.info;
    end
end

methods (Access=protected)
    buildproblem(obj);
end

end
