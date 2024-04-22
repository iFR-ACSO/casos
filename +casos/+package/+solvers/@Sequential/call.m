function argout = call(obj,argin)
    % Call sequential sum-of-squares
    
    import casos.package.UnifiedReturnStatus
    
    % get arguments
    args = argin;
    
    % initialize iteration info
    info = cell(1,obj.opts.max_iter);
    
    % first initial guess must be provided by user!
    xi_k = args{1};

    % check if we have parameter from the original nonlinear problem
    if ~is_zero(args{2})
        p0   = args{2};
    else
        p0 = [];
    end
    
    if obj.opts.verbose 
        disp('Start sequential SOS...')
        fprintf('%-10s%-15s%-15s%-15s\n', 'Iteration', 'abs(f)', 'abs(primal)', 'abs(dual)');
        fprintf('-------------------------------------------------\n');
    end

    % solve nonconvex SOS problem via sequence of convex SOS problem
    for i = 1:obj.opts.max_iter
        
        % initial guesss; actually not used but just for completness
        args{1} = xi_k;
    
        % set parameter to convex problem
        args{2}  = [p0; xi_k];
    
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
                
                constraint_vio = sqrt(full(casadi.DM( pnorm2( obj.constraintFun(xi_plus,p0) - sol{3} ) )));
                                
                % double(obj.cost_fun(xi_plus))
                cost      = double(sol{2});
                
  
            otherwise  
    
                if i > 1
    
                   % store iteration info
                   info(i+1:end) = [];
                   obj.info.iter = info;
                    
                   % take solution from previous iteration
                   argout = sol_old;
    
                  return
    
                else
    
                    % error: failed
                    obj.status = UnifiedReturnStatus.SOLVER_RET_NAN;
                    assert(~obj.opts.error_on_fail,'Convex optimization run into numerical errors.')
    
                end
      
            
        end
    

            % filter
            if i == 1
                Filter = [cost , constraint_vio];
                
            else
                % check new iterate
                % if Filter(:,2) < constraint_vio && Filter(:,1) < cost
                    % if both new values are larger then reject
                    % filterReject = 1;
                % else
                    % accept to filter
                    % filterReject = 0;
                    Filter(i,1) = cost;
                    Filter(i,2) = constraint_vio;
                % end
            end
    
            % check convergence
            if i > 1 
                    % select an algorithm to perform the linesearch
                    switch(obj.opts.line_search) 
                        case 'fminbnd'
                                dopt = line_search_fminbnd(obj, ...
                                                           p0,...
                                                           xi_k, ...
                                                           xi_plus, ...
                                                           dual_plus );
                        case 'bisection'
                                dopt = bisection_minimization(obj, ...
                                                              xi_k, ...
                                                              xi_plus, ...
                                                              dual_plus );
                        case 'polySolver'
                                sol = obj.lineSearch('p',[xi_k; xi_plus; dual_plus],...
                                                     'lbx',0.1,...
                                                     'ubx',1);

                                dopt = double(sol.x);


                        case 'none'
                            dopt = 1;
                    end

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
                sol_old = sol;

                else
    
                 % set current solution as previous solution for next iteration
                 xi_k = xi_plus;
                 dual_k = dual_plus;
                 sol_old = sol;
            end % end-if
    
    
    end % end for-loop sequential sos
end % end of function
