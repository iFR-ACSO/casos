classdef (Sealed) MosekInterface < casos.package.solvers.ConicSolver
% Interface for conic solver MOSEK.

properties (Access=protected)
    fhan;
    ghan;
    barv;
    cone;
end

properties (Access=private)
    info = struct;
end

properties (Constant, Access=protected)
    mosek_options = [casos.package.solvers.ConicSolver.conic_options
        {'mosek_param', 'Parameters to be passed to MOSEK.'
         'mosek_echo',  'Verbosity level passed to MOSEK (default: 0).'}
    ];
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.MosekInterface.mosek_options;
    end
end

methods
%     init(obj);
    argout = eval(obj,argin);

    function obj = MosekInterface(name,conic,varargin)
        % Construct MOSEK interface.
        obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});

        % default options
        if ~isfield(obj.opts,'mosek_param'), obj.opts.mosek_param = struct; end
        if ~isfield(obj.opts,'mosek_echo'), obj.opts.mosek_echo = 0; end
    end

    function s = stats(obj)
        % Return stats.
        s = obj.info;
        s = addfield(obj.status,s);
    end

    function sp = get_sparsity_in(obj,i)
        % Return sparsity pattern.
        if (i == 0)
            % Hessian pattern
            sp = sparsity(obj.args_in.h);

        else
            sp = get_sparsity_in@casos.package.solvers.ConicSolver(obj,i);
        end
    end
end

methods (Access=protected)
    buildproblem(obj);
end

end
