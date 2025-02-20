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

 
        % options to SCS
        opts = obj.opts.Clarabel;

        
        % Initialize an empty cell array for cones
        cones = cell(1, 1); 

        % Add zero cones for the Clarabel wrapper
        if K.z > 0
            cones{end+1} = zeroConeT(K.z);
        end

        % Add second order cones
        if K.q > 0
            cones{end+1} = SecondOrderConeT(K.q);
        end

        % Stack PSD cones
        num_psd = length(K.s);
        if num_psd > 0
            psd_cones = arrayfun(@PSDTriangleConeT, K.s, 'UniformOutput', false);
            % concatenate
            cones = [cones, psd_cones']; 
        end

        % Convert cell array to a proper structure 
        cones = vertcat(cones{:});
        
        % call clarabel mex
        sol = clarabel_mex(data.P,data.c,data.A,data.b,cones,opts);
        
        % extract primal, dual and slack variables
        x  = sol.x;
        z = sol.z;
        s = sol.s;
        
        % extract infos
        clarabel_info.solve_time = sol.solve_time;
        clarabel_info.status     = sol.status;
        clarabel_info.iterations = sol.iterations;
        clarabel_info.r_prim     = sol.r_prim;
        clarabel_info.r_dual     = sol.r_dual;
       
        obj.info.clarabel_info = clarabel_info;
       
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
