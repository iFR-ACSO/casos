classdef (Sealed) CDCSInterface < casos.package.solvers.ConicSolver
% Interface for conic solver Cdcs.

properties (Access=protected)
    fhan;
    ghan;
    cone;
end

properties (Access=private)
    info = struct;
end

properties (Constant, Access=protected)
    Cdcs_options = [casos.package.solvers.ConicSolver.conic_options
        {'Cdcs', 'Options to be passed to Cdcs.'}
    ];
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.CDCSInterface.Cdcs_options;
    end
end

methods
%     init(obj);
    argout = eval(obj,argin);

    function obj = CDCSInterface(name,conic,varargin)
        % Construct Cdcs interface.
        obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});

        % default options
        if ~isfield(obj.opts,'cdcs'), obj.opts.cdcs = []; end
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
