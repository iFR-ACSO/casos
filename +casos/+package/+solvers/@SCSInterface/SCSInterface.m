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

properties (Constant, Access=protected)
    scs_options = [casos.package.solvers.ConicSolver.conic_options
        {'scs', 'Options to be passed to SCS.'}
    ];
end

methods (Static, Access=protected)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SCSInterface.scs_options;
    end
end
methods
%     init(obj);
    argout = eval(obj,argin);

    function obj = SCSInterface(name,conic,varargin)
        % Construct SCS interface.
        obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});

        % default options
        if ~isfield(obj.opts,'scs'), obj.opts.scs = []; end
    end

    function s = stats(obj)
        % Return stats.
        s = obj.info;
        s = addfield(obj.status,s);
    end
end

methods (Access=protected)
    buildproblem(obj);
end

end
