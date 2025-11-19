classdef (Sealed) SCSInterface < casos.package.solvers.AbstractSCSInterface
% Interface for conic solver SCS.

properties (Constant, Access=protected)
    scs_options = [casos.package.solvers.ConicSolver.conic_options
        {'scs', 'Options to be passed to SCS.'}
    ];
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SCSInterface.scs_options;
    end
end

methods
    function obj = SCSInterface(name,conic,varargin)
        % Construct SCS interface.
        obj@casos.package.solvers.AbstractSCSInterface(name,conic,varargin{:});

        % default options
        if ~isfield(obj.opts,'scs'), obj.opts.scs = []; end
    end
end

methods (Static, Access=protected)
    %% Static helper functions
    function S = sparsity_triangular(k)
        % Return lower-triangular sparsity pattern.
        S = casadi.Sparsity.lower(k);
    end
end

methods (Access=protected)
    function [x,y,s] = call_solver(obj,data,K)
        % Call SCS solver.

        % options to SCS
        opts = obj.opts.scs;
        % disable verbosity by default
        if ~isfield(opts,'verbose'), opts.verbose = 0; end

        % call SCS
        [x,y,s,obj.solver_stats] = scs(data,K,opts);

        % check solver status
        status_val = obj.solver_stats.status_val;
        if status_val == -1
            % primal unbounded / dual infeasible
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
            assert(~obj.opts.error_on_fail,'Conic problem is dual infeasible.')
        elseif status_val == -2
            % primal infeasible / dual unbounded
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
            assert(~obj.opts.error_on_fail,'Conic problem is primal infeasible.')
        elseif ismember(status_val, [2 -6 -7])
            % inaccurate solution
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_NAN;
            assert(~obj.opts.error_on_fail,'Optimizer did not reach desired accuracy (Status: %s).', obj.solver_stats.status)
        elseif status_val < -2
            % failure
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
            assert(~obj.opts.error_on_fail,'Optimizer failed (Status: %s).', obj.solver_stats.status)
        else
            % success
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
        end
    end
end

end
