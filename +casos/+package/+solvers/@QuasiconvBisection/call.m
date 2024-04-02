function argout = call(obj,argin)
% Call bisection to solve quasiconvex SOS problem.

import casos.package.UnifiedReturnStatus

args = argin;

% prepare quasiconvex parameter 
if isscalar(argin{2})
    qcpar = repmat(argin{2},get_size_in(obj,1));
else
    qcpar = argin{2};
end

% substitute quasiconvex parameter
tvar = casos.PS.sym('t');
sossolver = substitute(obj.sossolver,'p',tvar,[qcpar; tvar]);

% initialize confidence intervals
interval = obj.qc_sign*obj.opts.conf_interval;

% initialize iteration info
info = cell(1,obj.opts.max_iter);

% last feasible solution
feas_sol = [];

printf(obj.log,'debug','====== Start Bisection =======\n\n');

for i=1:length(info)  
    % bisection
    ttry = mean(interval);

    % set parameter to convex problem
    args{2} = ttry;

    % evaluate convex SOS problem
    sol = call(sossolver, args);

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
            feas = 1;

        case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
            % infeasible: reset lower bound
            interval(1) = ttry;
            feas = 0;

        otherwise
            % error: failed
            feas = 0;
            obj.status = UnifiedReturnStatus.SOLVER_RET_NAN;
            printf(obj.log,'error','Convex optimization run into numerical errors.');
            assert(~obj.opts.error_on_fail,'Convex optimization run into numerical errors.')
            interval(1) = ttry;
    end

    % display bisection iterations
    tchar = 'gamma_';
    printf(obj.log,'info',['iteration:  %d/%d \t' tchar 'lb  = %-4.4f \t ' tchar 'try = %-4.4f \t'],...
                            i,length(info),obj.qc_sign*round(interval(1),4),obj.qc_sign*round(ttry,4));
    printf(obj.log,'info',[tchar 'ub = %-4.4f \t Feas = %d \n'],obj.qc_sign*round(interval(2),4),feas);
   
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
            printf(obj.log,'error','Quasiconvex problems appears infeasible.');
            assert(~obj.opts.error_on_fail,'Quasiconvex problems appears infeasible.')
            % return last solution
            argout = sol;
        end

        % store iteration info
        info(i+1:end) = [];
        obj.info.iter = info;

        printf(obj.log,'debug','\n====== Finished Bisection ======\n'); 

        % terminate
        return
    end
end

% no convergence
obj.status = UnifiedReturnStatus.SOLVER_RET_LIMITED;
printf(obj.log,'error','Bisection exceeded maximum number of iterations.');
assert(~obj.opts.error_on_fail,'Bisection exceeded maximum number of iterations.')

% return last solution
argout = sol;

% store iteration info
obj.info.iter = info;

end
