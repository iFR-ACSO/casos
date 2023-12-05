function argout = call(obj,argin)
% Call bisection to solve quasiconvex SOS problem.

import casos.package.UnifiedReturnStatus

args = argin;

% initialize confidence intervals
interval = obj.qc_sign*obj.opts.conf_interval;

% initialize iteration info
info = cell(1,obj.opts.max_iter);

% last feasible solution
feas_sol = [];

for i=1:length(info)
    % bisection
    ttry = mean(interval);

    % set parameter to convex problem
    args{2} = [argin{2}; ttry];

    % evaluate convex SOS problem
    sol = call(obj.sossolver, args);

    % store iteration info
    info{i} = obj.sossolver.stats;

    % set value
    sol{2} = casos.PS(obj.qc_sign*ttry);

    % update confidence interval
    switch (info{i}.UNIFIED_RETURN_STATUS)
        case UnifiedReturnStatus.SOLVER_RET_SUCCESS
            % feasible: reset upper bound
            interval(2) = ttry;
            % store feasible solution
            feas_sol = sol;

        case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
            % infeasible: reset lower bound
            interval(1) = ttry;

        otherwise
            % error: failed
            obj.status = UnifiedReturnStatus.SOLVER_RET_NAN;
            assert(~obj.opts.error_on_fail,'Convex optimization run into numerical errors.')
    end

    % check convergence
    if abs(diff(interval)) <= obj.opts.tolerance_abs ...
        || abs(diff(interval)) <= obj.opts.tolerance_rel*norm(interval)
        % either absolute or relative tolerance satisfied
        if ~isempty(feas_sol)
            % some feasible solution found
            obj.status = UnifiedReturnStatus.SOLVER_RET_SUCCESS;
            % return last feasible solution
            argout = feas_sol;

        else
            % no feasible solution found
            obj.status = UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
            assert(~obj.opts.error_on_fail,'Quasiconvex problems appears infeasible.')
            % return last solution
            argout = sol;
        end

        % store iteration info
        info(i+1:end) = [];
        obj.info.iter = info;

        % terminate
        return
    end
end

% no convergence
obj.status = UnifiedReturnStatus.SOLVER_RET_LIMITED;
assert(~obj.opts.error_on_fail,'Bisection exceeded maximum number of iterations.')

% return last solution
argout = sol;

% store iteration info
obj.info.iter = info;

end
