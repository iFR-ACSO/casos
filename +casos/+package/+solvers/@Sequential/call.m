function argout = call(obj,argin)
% Call sequential sum-of-squares

import casos.package.UnifiedReturnStatus

% get arguments
args = argin;

% initialize iteration info
info = cell(1,obj.opts.max_iter);

% first initial guess must be provided by user!
xk = args{1};

% do bisection
for i = 1:obj.opts.max_iter
    
    % initial guesss
    args{1} = xk;

    % set parameter to convex problem
    args{2}  = xk;

    % evaluate convex SOS problem
    sol = call(obj.sossolver, args);

    % store iteration info
    info{i} = obj.sossolver.stats;


    % check solution
    switch (info{i}.UNIFIED_RETURN_STATUS)
        case UnifiedReturnStatus.SOLVER_RET_SUCCESS
            
            % set initial guess and finally parameter; 
            % sos constraint reads  delta x = f(x0) + df/dx(x0)*(x-x0);
            % where delta x = xk+1 - xk hence xk+1 = delta_x + xk 
            xk1 = sol{1} + xk;
            sol_old = sol;
            % execute line search
            solLS = obj.lineSearch('p',[xk;xk1;sol{4};sol{5}], ...
                                   'lbx',0, ...
                                   'ubx',1);
            dopt = solLS.x;
         
        otherwise  
        % case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
            
            % Comments: How do we deal with solver infeasibility?
            if i > 1

               % store iteration info
               info(i+1:end) = [];
               obj.info.iter = info;

               % we take the last feasible solution
               xk1 = sol_old{1} + xk;
                
               sol_old{1} = xk1;

               argout = sol_old;

              return

            else

            % error: failed
            obj.status = UnifiedReturnStatus.SOLVER_RET_NAN;
            assert(~obj.opts.error_on_fail,'Convex optimization run into numerical errors.')
            end
        % otherwise
            % % error: failed
            % obj.status = UnifiedReturnStatus.SOLVER_RET_NAN;
            % assert(~obj.opts.error_on_fail,'Convex optimization run into numerical errors.')
        
    end


   % store iteration info
   info(i+1:end) = [];
   obj.info.iter = info;

    % check convergence
    if i > 1 
        
        % currently only primal variables
        xk1 = dopt*xk1 + (1-dopt)*xk;

        if full( casadi.DM(pnorm2(xk1-xk)) ) < obj.opts.tolerance_abs

        % adjust last solution similar to iteration i.e. overwrite
        % optimization solution 
        sol{1} = xk1;

        argout = sol;

        % terminate
        return

        end % end if


    end % end if

 % set current solution as previous solution for next iteration
 xk = xk1;

end % end for-loop

end % end of function
