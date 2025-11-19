classdef (Sealed) SedumiInterface < casos.package.solvers.ConicSolver
% Interface for conic solver SeDuMi.

properties (Access=protected)
    fhan;
    ghan;
    cone;

    solver_info  = struct;
    solver_stats = struct;
end

properties (Constant, Access=protected)
    sedumi_options = [casos.package.solvers.ConicSolver.conic_options
        {'sedumi', 'Options to be passed to SeDuMi.'}
    ];
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SedumiInterface.sedumi_options;
    end
end

methods
%     init(obj);
    argout = eval(obj,argin);

    function obj = SedumiInterface(name,conic,varargin)
        % Construct SeDuMi interface.
        obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});

        % default options
        if ~isfield(obj.opts,'sedumi'), obj.opts.sedumi = []; end
    end
    
    function s = stats(obj)
        % Return stats.
        s = obj.solver_stats;
        s = addfield(obj.status,s);
    end
end

methods (Access=protected)
    buildproblem(obj);
end

end
