classdef (Sealed) SedumiInterface < casos.package.solvers.ConicSolver
% Interface for conic solver SeDuMi.

properties (Access=protected)
    fhan;
    ghan;
    cone;
end

methods
%     init(obj);
    argout = eval(obj,argin);

    function obj = SedumiInterface(name,conic,varargin)
        % Construct SeDuMi interface.
        obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});
    end
end

methods (Access=protected)
    buildproblem(obj);
end

end
