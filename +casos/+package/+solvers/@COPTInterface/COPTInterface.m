classdef (Sealed) COPTInterface < casos.package.solvers.ConicSolver
% Interface for conic solver COPT.

properties (Access=protected)
    fhan;
    ghan;
    cone;
end

properties (Access=private)
    info = struct;
end

properties (Constant, Access=protected)
    copt_options = [casos.package.solvers.ConicSolver.conic_options
        {'copt', 'Options to be passed to COPT.'}
    ];
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.COPTInterface.copt_options;
    end
end

methods
%     init(obj);
    argout = eval(obj,argin);

    function obj = COPTInterface(name,conic,varargin)
        % Construct COPT interface.
        obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});

        % default options
        if ~isfield(obj.opts,'copt'), obj.opts.copt = []; end
    end

    function s = stats(obj)
        % Return stats.
        s = obj.info;
        s = addfield(obj.status,s);
    end
end

methods (Access=protected)
    buildproblem(obj);
    
    % vectorization of the semidefinite cone
    [v,i,j,k,l] = sdp_vec_upper(obj, M,Ks,scale,dim)
end

end

