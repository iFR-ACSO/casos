function argout = call(obj,argin)
    % Call sequential sum-of-squares
    
    import casos.package.UnifiedReturnStatus
    
    % get arguments
    args = argin;
    
    % initialize iteration info
    info = cell(1,obj.opts.max_iter);
    
    % first initial guess must be provided by user!
    xi_k = args{1};
    
    if obj.opts.verbose 
    disp('Start sequential SOS...')
    fprintf('%-10s%-15s%-15s%-15s\n', 'Iteration', 'abs(f)', 'abs(primal)', 'abs(dual)');
    fprintf('-------------------------------------------------\n');
    end

    % solve nonconvex SOS problem via sequence of convex SOS problem
    for i = 1:obj.opts.max_iter
        
        % initial guesss
        args{1} = xi_k;
    
        % set parameter to convex problem
        args{2}  = xi_k;
    
        % evaluate convex SOS problem
        sol = call(obj.sossolver, args);
    
        % store iteration info
        info{i} = obj.sossolver.stats;
    
    
        % check solution
        switch (info{i}.UNIFIED_RETURN_STATUS)
            case UnifiedReturnStatus.SOLVER_RET_SUCCESS
                
                % set initial guess and finally parameter; 
                % sos constraint reads  delta x = f(x0) + df/dx(x0)*(x-x0);
                % where delta x = xi_k+1 - xi_k hence xi_k+1 = delta_x + xi_k 
                xi_plus   = sol{1} + xi_k;
    
                dual_plus = sol{5};
    
                cost = double(sol{2});
                
                % store old solution
                sol_old    = sol;
                sol_old{1} =  xi_plus;
    
            otherwise  
    
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
    
    
            % check convergence
            if i > 1 
                    
                    dopt = line_search_fminbnd(obj,  xi_k,xi_plus,sol);
     
                    % update primal and dual solution
                    xi_k1    = dopt*xi_plus    + (1-dopt)*xi_k;
                    duals_k1 = dopt*dual_plus  + (1-dopt)*dual_k;
                     
    
                 if obj.opts.verbose 
                    fprintf('%-10d%-15.4f%-15.4f%-15.4f\n',...
                            i, cost ,...
                            full(casadi.DM(pnorm2(xi_k1 - xi_k))), full(casadi.DM(pnorm2(duals_k1 - dual_k))));
                end
    
                % check convergence or if maximum number of iterations are reached    
                if i == obj.opts.max_iter || ...
                   full( casadi.DM( pnorm2(xi_k1 - xi_k)) )       < obj.opts.tolerance_abs && ...
                   full( casadi.DM( pnorm2(duals_k1 - dual_k) ) ) < obj.opts.tolerance_rel*full( casadi.DM( pnorm2( duals_k1 ) ) ) 
                  
                   % store iteration info
                   info(i+1:end) = [];
                   obj.info.iter = info;
            
                    % adjust last solution similar to iteration i.e. overwrite
                    % optimization solution 
                    sol{1} = xi_k1;
                    sol{5} = duals_k1;
            
                    argout = sol;
            
                    % terminate
                    return
        
                end % end-if 

                % set current solution as previous solution for next iteration
                xi_k      = xi_k1;
                dual_k    = duals_k1;
         
            else
    
                 % set current solution as previous solution for next iteration
                 xi_k = xi_plus;
                 dual_k = dual_plus;
    
            end % end-if
    
    
    end % end for-loop sequential sos
end % end of function
