classdef (Sealed) ClarabelInterface < casos.package.solvers.ConicSolver
% Interface for conic solver Clarabel.

properties (Access=protected)
    fhan;
    ghan;
    cone;

    cones;
end

properties (Access=private)
    info;
end

properties (Constant, Access=protected)
    Clarabel_options = [casos.package.solvers.ConicSolver.conic_options
        {'Clarabel', 'Options to be passed to Clarabel.'}
    ];
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.ClarabelInterface.Clarabel_options;
    end
end
methods
%     init(obj);
    argout = eval(obj,argin);

    function obj = ClarabelInterface(name,conic,varargin)
        % Construct Clarabel interface.
        obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});

        % default options; see default setting structure of clarabel matlab
        % interface
        if ~isfield(obj.opts,'Clarabel'), obj.opts.Clarabel = DefaultSettings; end
    end

    function s = stats(obj)
        % Return stats.
        s = obj.info;
        s = addfield(obj.status,s);
    end
end

methods (Static, Access=protected)
    %% Static helper functions
    % MOSEK/Clarabel-style for the semidefinite cone
    V = sdp_vec(M,varargin);
    M = sdp_mat(V,varargin);
end

methods (Access=protected)
    buildproblem(obj);
end

end
