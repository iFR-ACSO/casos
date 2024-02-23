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
            
           %  xCoeff_sol     = poly2basis(xk1);
           %  lamCoeff_sol   = poly2basis(sol{5});
           % % 
           % nabla_x = full(casadi.DM(obj.nabla_x_fun(xk, sol{1}, sol{5})))
           % nabla_l = full(casadi.DM(obj.nabla_lam_fun(xk, sol{1}, sol{5})))

            % store old solution
            sol_old    = sol;
            sol_old{1} =  xk1;

        otherwise  
        % case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
            
            % Comments: How do we deal with solver infeasibility?
            if i > 1

               % store iteration info
               info(i+1:end) = [];
               obj.info.iter = info;

               argout = sol_old;

              return

            else

            % error: failed
            obj.status = UnifiedReturnStatus.SOLVER_RET_NAN;
            assert(~obj.opts.error_on_fail,'Convex optimization run into numerical errors.')

            end
  
        
    end


       % plot merit function
       % figure(1)
       % % clf
       % plot(linspace(0,1,10),arrayfun(@(d) double(obj.Merit(xk.*(1-d)+xk1.*d,sol{5})),...
       %     linspace(0,1,10)))

        

        
        % check convergence
        if i > 1 
                
                dopt = bisection_minimization(obj, xk, xk1, sol);

           
                % solLS = obj.lineSearch('p',[xk;xk1;sol{5}], ...
                %                        'lbx',0, ...
                %                        'ubx',1);


                % dopt = solLS.x;
        
                %update primal
                xk1   = dopt*xk1 + (1-dopt)*xk;
                duals =  dopt*sol{5} + (1-dopt)*sol{5};

            if full( casadi.DM(pnorm2(xk1-xk)) ) < obj.opts.tolerance_abs %&& ...
               full( casadi.DM(pnorm2(duals-sol{5})) ) < obj.opts.tolerance_rel*full( casadi.DM(pnorm2(duals))) 
               % store iteration info
               info(i+1:end) = [];
               obj.info.iter = info;
        
                % adjust last solution similar to iteration i.e. overwrite
                % optimization solution 
                sol{1} = xk1;
        
                argout = sol;
        
                % terminate
                return
    
            end % end if convergence check

            xk = xk1;
        else
             % set current solution as previous solution for next iteration
             xk = xk1;
             % dk = sol{5};
        end


end % end for-loop sequential sos

end % end of function
