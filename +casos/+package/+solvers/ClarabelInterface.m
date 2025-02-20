classdef (Sealed) ClarabelInterface < casos.package.solvers.AbstractSCSInterface
% Interface for conic solver Clarabel.

properties (Constant, Access=protected)
    clarabel_options = [casos.package.solvers.ConicSolver.conic_options
        {'clarabel', 'Options to be passed to Clarabel.'}
    ];
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.ClarabelInterface.clarabel_options;
    end
end

methods
    function obj = ClarabelInterface(name,conic,varargin)
        % Construct Clarabel interface.
        obj@casos.package.solvers.AbstractSCSInterface(name,conic,varargin{:});

        % default options
        % see default setting structure of Clarabel MATLAB interface
        if ~isfield(obj.opts,'Clarabel'), obj.opts.Clarabel = DefaultSettings; end
    end
end

methods (Static, Access=protected)
    %% Static helper functions
    function S = sparsity_triangular(k)
        % Return upper-triangular sparsity pattern.
        S = casadi.Sparsity.upper(k);
    end
end

methods (Access=protected)
    function [x,z,s] = call_solver(obj,data,K)
        % Call Clarabel solver.

        cones = [];

        % options to SCS
        opts = obj.opts.clarabel;

        % add zero cones
        if K.z > 0
            cones(end+1,:) = zeroConeT(K.z);
        end
        
        % add second order cones
        if K.q > 0
            cones(end+1,:) = SecondOrderConeT(K.q);
        end
        
        % stack PSD cones
        if ~isempty(K.s)
            psd_cones = arrayfun(@PSDTriangleConeT, K.s, 'UniformOutput', false);
            % concatenate
            cones = vertcat(cones, psd_cones{:}); 
        end
        
        % call clarabel mex
        sol = clarabel_mex(data.P,data.c,data.A,data.b,cones,opts);
        
        % extract primal, dual and slack variables
        x  = sol.x;
        z = sol.z;
        s = sol.s;

        % get solution status
        if strcmp(sol.status,'Solved')
            % success
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
        elseif strcmp(sol.status,'AlmostSolved')
            % success because we are feasible but did not solve until thresholds
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
        elseif strcmp(sol.status,'Unsolved')
            % unsolved
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        elseif strcmp(sol.status,'MaxIterations')
            % max iterations reached
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        elseif strcmp(sol.status,'MaxTime')
            % max time reached
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        elseif strcmp(sol.status,'InsufficientProgress')
            % InsufficientProgress
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        elseif strcmp(sol.status,'Unknown')
            % InsufficientProgress
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
        elseif strcmp(sol.status,'NumericalError')
            % InsufficientProgress
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        else
            % primal unbounded / dual infeasible
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
            assert(~obj.opts.error_on_fail,['Conic problem is ' sol.status])
        end
    end
end

end
